# Sander-Python API

The sander API is exposed to Python via the `pysander` bindings implemented in this directory.
The work here leverages the sander C-API (itself built on the Fortran API) that is compiled as part of the sander build.
Note that there are two components here:

1. A Python wrapper (found in `sander/__init__.py`) that does basic type and value checking of input parameters.
1. A C-API wrapper (found in `sander/src`) that contains the Python type (struct) definitions (`pysandermoduletypes.c`) and wrappers around the C-API (`pysandermodule.c`).

## Adding a new type

All of the types are basic C `struct`-like types that mirror the relevant data structures in the C-API.
These are defined in `pysandermoduletypes.c`, and the existing types should serve as example templates in the case additional types are needed in the future.

Here are the steps needed:

### Define the types in `pysandermoduletypes.c`

A `struct`-like type requires 5 separate definitions:

#### Define the data structure

Each member of the data structure should be of type `PyObject*`:

```C
typedef struct {
    PyObject_HEAD
    PyObject *arg1;     // int
    PyObject *arg2;     // float
    PyObject *arg3;     // str
} pysander_FooOptions;
```

#### Define the allocator

When allocating the data structure, we must make sure to properly handle the reference counts to avoid leaking or segfaulting.

```C
static PyObject *
pysander_FooOptions_new(PyTypeObject *type) {
    pysander_FooOptions *self;
    self = (pysander_FooOptions *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->arg1 = PyLong_FromLong(0);
        self->arg2 = PyFloat_FromDouble(0.0);

        // There is a helpful ASSIGN_STRING macro for dealing with unicode
        ASSIGN_STRING(arg3, "default");
    }

    return (PyObject *) self;
}
```

#### Define the member structure

This is useful for describing an attribute that corresponds to a C struct member.
You can see more in the [official Python documentation.](https://docs.python.org/3/c-api/structures.html#c.PyMemberDef)

```C
static PyMemberDef pysander_FooOptionMembers[] = {
    {"arg1", T_OBJECT_EX, offsetof(pysander_FooOptions, arg1), 0, "First argument description"},
    {"arg2", T_OBJECT_EX, offsetof(pysander_FooOptions, arg2), 0, "Second argument description"},
    {"arg3", T_OBJECT_EX, offsetof(pysander_FooOptions, arg3), 0, "Third argument description"},

    {NULL}
};
```

#### Define the deallocator

It is necessary to decrement the refcounter on each variable that is freed or else memory will be leaked.
This looks like the following:

```C
static void
pysander_FooOptions_dealloc(pysander_FooOptions* self) {
    Py_DECREF(self->arg1);
    Py_DECREF(self->arg2);
    Py_DECREF(self->arg3);

    Py_TYPE(self)->tp_free((PyObject *)self);
}
```

#### Define the type itself

This is the type definition available in Python itself. It looks like:

```C
static PyTypeObject pysander_FooOptionsType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "sander.pysander.FooOptions",   // tp_name
    sizeof(pysander_FooOptions),    // tp_basicsize
    0,                              // tp_itemsize
    (destructor)pysander_FooOptions_dealloc, // tp_dealloc
    0,                              // tp_print
    0,                              // tp_getattr
    0,                              // tp_setattr
    0,                              // tp_compare
    0,                              // tp_repr
    0,                              // tp_as_number
    0,                              // tp_as_sequence
    0,                              // tp_as_mapping
    0,                              // tp_hash
    0,                              // tp_call
    0,                              // tp_str
    0,                              // tp_getattro
    0,                              // tp_setattro
    0,                              // tp_as_buffer
    Py_TPFLAGS_DEFAULT,             // tp_flags
    "Description of Foo options",   // tp_doc
    0,		                        // tp_traverse
    0,		                        // tp_clear
    0,		                        // tp_richcompare
    0,		                        // tp_weaklistoffset
    0,		                        // tp_iter
    0,		                        // tp_iternext
    0,                              // tp_methods
    pysander_FooOptionMembers,      // tp_members
    0,                              // tp_getset
    0,                              // tp_base
    0,                              // tp_dict
    0,                              // tp_descr_get
    0,                              // tp_descr_set
    0,                              // tp_dictoffset
    0,                              // tp_init
    0,                              // tp_alloc
    (newfunc)pysander_FooOptions_new,// tp_new
};
```

Notice the reference to the functions we defined before. As a result, this part must be defined *after* the others.

### Register the new type in the extension module

Once you have defined your new data structure above, you have to make it available in the namespace of the extension module.
To do this, the type must be registered with the module.
This is done for all types at the end of the `PyInit_pysander` function in `pysandermodule.c`:

```C
PyMODINIT_FUNC
PyInit_pysander(void) {
    // Type declarations
    if (PyType_Ready(&pysander_InputOptionsType) < 0)
        return NULL;
    if (PyType_Ready(&pysander_EnergyTermsType) < 0)
        return NULL;
    if (PyType_Ready(&pysander_QmInputOptionsType) < 0)
        return NULL;

    // All ready checks should go here
    if (PyType_Ready(&pysander_FooOptionsType) < 0)
        return NULL;

    PyObject* m = PyModule_Create(&moduledef);

    // Now add the types
    Py_INCREF(&pysander_InputOptionsType);
    PyModule_AddObject(m, "InputOptions", (PyObject*) &pysander_InputOptionsType);
    Py_INCREF(&pysander_EnergyTermsType);
    PyModule_AddObject(m, "EnergyTerms", (PyObject *) &pysander_EnergyTermsType);
    Py_INCREF(&pysander_QmInputOptionsType);
    PyModule_AddObject(m, "QmInputOptions", (PyObject *) &pysander_QmInputOptionsType);

    // INCREF the type definition and then add it to the module namespace
    Py_INCREF(&pysander_FooOptionsType);
    PyModule_AddObject(m, "FooOptions", (PyObject *) &pysander_FooOptionsType);

    return m;
}
```

## Expanding the sander setup variables

### Expanding the extension module (`pysandermodule.c`)

If a new set of options needs to be supplied to the C setup routine (and the relevant changes have been made to the C-API of sander),
then you can pass in the new input types by expanding the argument processing code in `pysander_setup()` in `pysandermodule.c`:

```C
    if (!PyArg_ParseTuple(args, "sOOOO", &prmtop, &arg2, &arg3, &arg4, &arg5))
        return NULL;
```

You can find the documentation for `PyArg_ParseTuple` [here](https://docs.python.org/3/c-api/arg.html#c.PyArg_ParseTuple).
As a quick description, the current call parses 5 arguments: a `const char*` (`s`) and 4 `PyObject*`s (`O`).
Optional arguments (such as `arg5`, which is the optional QM arguments) should be set to `None` when not provided, and
can be compared to `Py_None` in the logic flow of `pysander_setup`.
The existing code provides examples of this code flow.

To provide the `FooOptions` here, the `PyArg_ParseTuple` call should be expanded to:

```C
    // Don't forget to declare arg6 above
    if (!PyArg_ParseTuple(args, "sOOOOO", &prmtop, &arg2, &arg3, &arg4, &arg5, &arg6))
        return NULL;
```

Then add the necessary code handle `arg6` if provided (or is set to `Py_None`) below.

### Adding support in the Python wrapper

The Python wrapper is defined in `__init__.py` (look for the call to `_pys.setup`).
The API users will call is defined *here*, and the interface (including defaulting to `None`) should be handled in `__init__.py`.
Also note how the extension module type definitions are exposed via `__init__.py` - you must also make sure your new `FooOptions` type is exposed here as well if users should have access to it.
