#!/bin/sh

. `dirname $0`/test_check.sh

date_string=`date +%Y-%m-%d_%H-%M-%S`
logdir="${AMBERHOME}/logs/test_at_openmp"
logprefix="${logdir}/${date_string}"
logfile="${logprefix}.log"
difffile="${logprefix}.diff"

dir=`dirname \`pwd\``
check_environment `dirname $dir` openmp

mkdir -p ${logdir}
(make -k -f Makefile test.openmp2 2>&1) | tee ${logfile}
(make -k -f Makefile finished 2>&1) | tee -a ${logfile}

passed_count=`grep PASS ${logfile} | wc -l`
questionable_count=`grep "FAILURE:" ${logfile} | wc -l`
questionable_ignored_count=`grep -e "FAILURE: (ignored)" ${logfile} | wc -l`
error_count=`grep "Program error" ${logfile} | grep -v "echo" | wc -l`

echo "${passed_count} file comparisons passed" | tee -a ${logfile}
if [ ${questionable_count} -eq 0 ]; then
    echo "${questionable_count} file comparisons failed" | tee -a ${logfile}
else
    echo "${questionable_count} file comparisons failed (${questionable_ignored_count}" \
         "of which can be ignored)" | tee -a ${logfile}
fi
echo "${error_count} tests experienced errors" | tee -a ${logfile}

echo "Test log file saved as ${logfile}" | tee -a ${logfile}

if [ -f TEST_FAILURES.diff ]; then
   mv TEST_FAILURES.diff ${difffile}
   if [ -f ${AMBERHOME}/test/TEST_FAILURES.diff ]; then
      cat ${AMBERHOME}/test/TEST_FAILURES.diff >> ${difffile}
      /bin/rm ${AMBERHOME}/test/TEST_FAILURES.diff
   fi
   echo "Test diffs file saved as ${difffile}" | tee -a ${logfile}
else
   if [ -f ${AMBERHOME}/test/TEST_FAILURES.diff ]; then
      mv ${AMBERHOME}/test/TEST_FAILURES.diff ${difffile}
      echo "Test diffs file saved as ${difffile}" | tee -a ${logfile}
   else
      echo "No test diffs to save!" | tee -a ${logfile}
   fi
fi

# save summary for later reporting:
tail -5 ${logfile} > ${logdir}/at_summary

if [ ${questionable_count} -gt ${questionable_ignored_count} -o ${error_count} -gt 0 ]; then
    # Make sure the exit code reflects the error count if we have comparison failures that
    # are not intended to be ignored
    exit 1
fi
