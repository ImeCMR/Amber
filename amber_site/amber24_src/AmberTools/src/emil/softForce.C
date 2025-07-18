#include "config.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


using namespace std;
#include "emil.h"



#define SF_FORCE_R   SF_MAX_R   //Force goes to zero at this distance
#define SF_FORCE_R2  SF_MAX_R2  //these are deffed in the main emil.h
 

//#define SOFT_FORCE_CONE
//#define SOFT_FORCE_SMOOTHCONE


//#define SF_WELL_H 180.0
//#define SF_EPS_MUL 1.35

//third order smoothstep function: derivatives are zero at 0,1.
//                                 values are 0,1       at 0,1.
inline double smoothStep_3( double x, double *dfdx ){
  double x2;
    x2  = x*x;
  *dfdx = 6*(x - x2);
  return ( -2*x2*x + 3*x2 );
}


//calculate energy only
double hsc::getSoftEnergyAll(Particle *p ) {
  
  Cell     *curCell;
  Particle *curPart;
  double    softE;
  
  softE = 0.0;
  
  curCell = p->cell;
  curPart = curCell->firstParticle;

  // Local cell
  while (curPart){    
      if( curPart != p )
        softE += getSoftEnergy( p, curPart );
    curPart = curPart->next;
  }

  
  // Neighbour cells
  for(int j=0; j<26; j++){
  
    curCell = p->cell->neighbours[j];
    curPart = curCell->firstParticle;
   
    while (curPart){
      softE += getSoftEnergy( p, curPart );
      curPart = curPart->next;
    }
  }
  
  return( softE );
}
  
double hsc::softForceAll(Particle *p, double scaleSoft, double *fij ) {

  Cell     *curCell;
  Particle *curPart;
  double    softE, delta;
  
  softE = 0.0;
  
  curCell = p->cell;
  curPart = curCell->firstParticle;

  // Local cell
  while (curPart){    
      if( curPart != p ){
        softE += softForce( p, curPart, fij, scaleSoft);
      }
    curPart = curPart->next;
  }

  
  // Neighbour cells
  for(int j=0; j<26; j++){
  
    curCell = p->cell->neighbours[j];
    curPart = curCell->firstParticle;
   
    while (curPart){
      softE  += softForce( p, curPart, fij, scaleSoft);
      curPart = curPart->next;
    }
  }
  
  return( softE );
      
} 

double hsc::getSoftEnergy(Particle *p1, Particle *p2 ) {

      double r2, rij[3]; 
      double r, E;

#ifdef PER_PAIR_SC
      int    pairIndex, pair;
      double myeps, ljA;

      //load the epsilon used in softcoring the interaction between this atom pair
      pairIndex = (LJ_atomTypes[p1->myId] - 1) * nLJ_atomTypes + LJ_atomTypes[p2->myId] - 1;
      pair      = LJ_pairTypes[pairIndex];
      myeps     = SF_EPS_MUL * scEpsilon( pair );
      ljA       = sc_ljA( pair );
#endif

      //find the distance
      if(periodic){
        for (int j=0; j<3; j++){
            rij[j] = p1->R[j] - p2->R[j];
            while (fabs(rij[j]) > box->halfx[j]) 
              rij[j] -= copysign(box->x[j], rij[j]);
        }
      }else{
        for (int j=0; j<3; j++){
            rij[j] = p1->R[j] - p2->R[j];
        } 
      }
      r2 = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
      
      //define "overlap" as approach closer than epsilon.
      if ( r2 >= SF_FORCE_R2 ) return 0.0;      

      if( r2 == 0.0 ){
          MERR << "Soft force for " << p1->myId << " and " << p2->myId << " r = " << sqrt(r2) << " E = " << SF_WELL_H << endl;
          return SF_WELL_H;
      }
          
      //if we are still here, there is a force
      r    = sqrt( r2 );
      //E     = SF_WELL_H * pow( (SF_FORCE_R - r)/SF_FORCE_R, SF_WELL_EXP );
      
#ifdef SOFT_FORCE_CONE
      E  = SF_WELL_H * (SF_FORCE_R - r)/SF_FORCE_R;
#elifdef SOFT_FORCE_SMOOTHCONE
      E  = 2 * SF_WELL_H * smoothStep_3( (SF_FORCE_R - 0.5*r)/SF_FORCE_R, &r2);
#else
      E  = SF_WELL_H * smoothStep_3( (SF_FORCE_R - r)/SF_FORCE_R, &r2);
#endif
      
      return( E );
}

double hsc::softForce(Particle *p1, Particle *p2, double *fij, double scaleSoft ) {

      double r2, rij[3]; 
      double r, E, f;

#ifdef PER_PAIR_SC
      int    pairIndex, pair;
      double myeps, ljA;

      //load the epsilon used in softcoring the interaction between this atom pair
      pairIndex = (LJ_atomTypes[p1->myId] - 1) * nLJ_atomTypes + LJ_atomTypes[p2->myId] - 1;
      pair      = LJ_pairTypes[pairIndex];
      myeps     = SF_EPS_MUL * scEpsilon( pair );
      ljA       = sc_ljA( pair );
#endif

      
      //find the distance
      if(periodic){
        for (int j=0; j<3; j++){
          rij[j] = p1->R[j] - p2->R[j];
          while (fabs(rij[j]) > box->halfx[j]) 
            rij[j] -= copysign(box->x[j], rij[j]);
          }
      }else{
        for (int j=0; j<3; j++){
          rij[j] = p1->R[j] - p2->R[j];
        } 
      }
      r2 = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
      
      //define "overlap" as approach closer than epsilon.
      if( r2 >= SF_FORCE_R2 ) return 0.0;
      
      if( r2 == 0.0 ){
          MERR   << "EMIL Warning, exact position overlap between particles " << p1->myId << " and " << p2->myId << endl; 
        *logFile << "esoft of " << scaleSoft * SF_WELL_H << " between "<< p1->myId << " and " << p2->myId << endl; 
          return SF_WELL_H;
      }
          
      //if we are still here, there is a force
      r    = sqrt( r2 );
          
      //repulsion is some power of R.
 //     E     = SF_WELL_H * pow( (SF_FORCE_R - r)/SF_FORCE_R, SF_WELL_EXP );
 //     f     = SF_WELL_EXP * E * (SF_FORCE_R / (SF_FORCE_R - r));//f is dE/dr

#ifdef SOFT_FORCE_CONE
      E   = SF_WELL_H * (SF_FORCE_R - r)/SF_FORCE_R;
      f   = SF_WELL_H / SF_FORCE_R;
#elifdef SOFT_FORCE_SMOOTHCONE
      E   = 2 * SF_WELL_H * smoothStep_3( (SF_FORCE_R - 0.5*r)/SF_FORCE_R, &f);
      f  *= 2 * SF_WELL_H;
      f  /= r;         
#else
      E   = SF_WELL_H * smoothStep_3( (SF_FORCE_R - r)/SF_FORCE_R, &f);
      f  *= SF_WELL_H;
      f  /= r;         
#endif

      //force on p1 due to p2.
      f    *= scaleSoft;    //scale force before adding, but do not scale energy yet. 
         
      fij[0] += f * rij[0]; //f is positive and vector rij points from p1 to p2.
      fij[1] += f * rij[1];
      fij[2] += f * rij[2];
       
      
      return( E );
}
      
bool hsc::checkOverlapSoft(Particle &p1, Particle &p2) {

      double rij[3], rijsq, epsilon;
      int    pairIndex, pair;
#if 0  
      //load the epsilon used in softcoring the interaction between this atom pair
      pairIndex = (LJ_atomTypes[p1.myId] - 1) * nLJ_atomTypes + LJ_atomTypes[p2.myId] - 1;
      pair      = LJ_pairTypes[pairIndex];
      epsilon   = SF_EPS_MUL * scEpsilon( pair );
      
      if(periodic){
        for (int j=0; j<3; j++){
          rij[j] = p1.R[j] - p2.R[j];
          while (fabs(rij[j]) > box->halfx[j]) 
            rij[j] -= copysign(box->x[j], rij[j]);
        }
      }else{
        for (int j=0; j<3; j++){
          rij[j] = p1.R[j] - p2.R[j];
        } 
      }
      rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
    
      
      //define "overlap" as approach closer than epsilon
      if (rijsq < epsilon * epsilon) return true;
      else 
#endif
        return true;
      
}
#undef SF_WELL_H
