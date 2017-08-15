#include<R_ext/Rdynload.h>
#ifndef R_R_H
#  include <R.h>
#endif
/*     subroutine advection(N, Y, dt, h, hint, v, Bcup, Bcdown,              &
     &       Yup, Ydown, VFint, VF, Aint, A, method, mode, split,            &
     &       dY, cu, it 
     
       subroutine advectvol(N, Y, dt, V, Vint, flow, Bcup, Bcdown,           &
     &                     Yup, Ydown,  method,mode,dY, cu, it)

      SUBROUTINE diffspher (Nr, Ntet, Nphi, C,                                   &
     &             Crup, Crdown, Ctetup, Ctetdown, Cphiup, Cphidown,             &
     &             Frup, Frdown, Ftetup, Ftetdown, Fphiup, Fphidown,             &
     &             Bcrup, Bcrdown, Bctetup, Bctetdown,                           &
     &             Bcphiup, Bcphidown, D_r, D_tet, D_phi,                        &
     &             r, tet, phi, FluxrUp, FluxrDown,                              &
     &             FluxtetUp, FluxtetDown, FluxphiUp, FluxPhiDown, dC) 

      SUBROUTINE diffcyl   (Nr, Ntet, Nz, C,                                     &
     &             Crup, Crdown, Ctetup, Ctetdown, Czup, Czdown,                 &
     &             Frup, Frdown, Ftetup, Ftetdown, Fzup, Fzdown,                 &
     &             Bcrup, Bcrdown, Bctetup, Bctetdown,                           &
     &             BCzUp, BCzDown, D_r, D_tet, D_z,                              &
     &             r, tet, z, FluxrUp, FluxrDown,                                &
     &             FluxtetUp, FluxtetDown, FluxzUp, FluxzDown, dC) 

     */

void F77_NAME(advection)(int*, double*, double*, double*, double*, double *, 
                         int*, int*, 
                         double*, double*, double*, double*, double*, double*, 
                         int*, int*, int*, double*, double*, int*);

void F77_NAME(advectvol)(int*, double*, double*, double*, double*, double *, 
                         int*, int*, 
                         double*, double*, int*, int*, double*, double*, int*);

void F77_NAME(diffcyl) (int*, int*, int*, double*,  
      double*, double*, double*, double*, double*, double*, 
      double*, double*, double*, double*, double*, double*, 
      int*, int*, int*, int*, int*, int*,  
      double*, double*, double*, double*, double*, double*,
      double*, double*, double*, double*, double*, double*, double*);

void F77_NAME(diffspher) (int*, int*, int*, double*,  
      double*, double*, double*, double*, double*, double*,
      double*, double*, double*, double*, double*, double*, 
      int*, int*, int*, int*, int*, int*,  
      double*, double*, double*, double*, double*, double*,
      double*, double*, double*, double*, double*, double*, double*);  

      
R_FortranMethodDef ReacfortranMethods[] = {
 {"advection",   (DL_FUNC) &F77_SUB(advection), 20},
 {"advectvol",   (DL_FUNC) &F77_SUB(advectvol), 15},
 {"diffcyl"  ,   (DL_FUNC) &F77_SUB(diffcyl),   35},
 {"diffspher" ,  (DL_FUNC) &F77_SUB(diffspher), 35},
 {NULL, NULL, 0}
};

void R_init_ReacTran(DllInfo *info) {
  R_registerRoutines (info, NULL, NULL, ReacfortranMethods, NULL);
  R_useDynamicSymbols(info, FALSE); // disable dynamic searching  
}
