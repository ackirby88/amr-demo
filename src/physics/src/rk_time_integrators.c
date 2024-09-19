/**
 * \file   rk_time_integrators.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Runge-Kutta ODE Integrators
 * 
 * Explicit Runge-Kutta methods available:
 *      * ERK1: one stage, first-order accuracy
 *      * ERK2: two stage, second-order accuracy
 *      * ERK3: three stage, third-order accuracy
 *      * ERK4: four stage, fourth-order accuracy
 *      * LSERK45: five stage, fourth-order accuracy, low storage
 * 
 * Created on August 23, 2018, 12:42 PM
 */

/* header files */
#include "rk_time_integrators.h"


void rk_1(  int* dim,
            int* unstructured,
            int* nfields,
            int* nelem_subgrid,
            int* ncell,
            int* nghost,
            int* nface,
            int* dof,
            int* cell_info,
            int* face_info,
            int* nsend,
            int* nrecv,
            int* send,
            int* recv,
            int* sizes,
            double* geom,
            double* Q,
            double* k1,
            double dt,
            MPI_Comm mpicomm_ext,
            double* sendbuf,
            MPI_Request* request,
            int* counts,
            double* resnorm,
            double* compute_time,
            double* residual_time){

    int i;
    double t1,t2;
    double beta1;
    double resnorm_in;
    
    
    //RK1
    beta1  = 1.0*dt;
    
    /* stage 1 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k1);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;

    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k1);
    for(i=0;i<*dof;++i){
        Q[i] = Q[i] + beta1*k1[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    MPI_Allreduce(&resnorm_in,resnorm,1,MPI_DOUBLE,MPI_SUM,mpicomm_ext);
    *resnorm = sqrt(*resnorm);
    
}


void rk_2(  int* dim,
            int* unstructured,
            int* nfields,
            int* nelem_subgrid,
            int* ncell,
            int* nghost,
            int* nface,
            int* dof,
            int* cell_info,
            int* face_info,
            int* nsend,
            int* nrecv,
            int* send,
            int* recv,
            int* sizes,
            double* geom,
            double* Q,
            double* Q_stage,
            double* k1,
            double* k2,
            double dt,
            MPI_Comm mpicomm_ext,
            double* sendbuf,
            MPI_Request* request,
            int* counts,
            double* resnorm,
            double* compute_time,
            double* residual_time){

    int i;
    double t1,t2;
    double alpha1;
    double beta1;
    double beta2;
    double resnorm_in;
    
    //RK2 SSP TVD
    alpha1 = 1.0*dt;
    beta1  = 0.5*dt;
    beta2  = 0.5*dt;
            
    
    //save Q
    for(i=0;i<*dof;++i){
        Q_stage[i] = Q[i];
    }
    
    
    /* stage 1 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k1);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k1);
    for(i=0;i<*dof;++i){
        Q[i] = Q_stage[i] + alpha1*k1[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    
    
    /* stage 2 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k2);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k2);
    for(i=0;i<*dof;++i){
        Q[i] = Q_stage[i] + beta1*k1[i] + beta2*k2[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    MPI_Allreduce(&resnorm_in,resnorm,1,MPI_DOUBLE,MPI_SUM,mpicomm_ext);
    *resnorm = sqrt(*resnorm);
    
}



void rk_3(  int* dim,
            int* unstructured,
            int* nfields,
            int* nelem_subgrid,
            int* ncell,
            int* nghost,
            int* nface,
            int* dof,
            int* cell_info,
            int* face_info,
            int* nsend,
            int* nrecv,
            int* send,
            int* recv,
            int* sizes,
            double* geom,
            double* Q,
            double* Q_stage,
            double* k1,
            double* k2,
            double* k3,
            double dt,
            MPI_Comm mpicomm_ext,
            double* sendbuf,
            MPI_Request* request,
            int* counts,
            double* resnorm,
            double* compute_time,
            double* residual_time){
 
    int i;
    double t1,t2;
    double alpha1;
    double alpha21;
    double alpha22;
    double beta1;
    double beta2;
    double beta3;
    double resnorm_in;
    
    
    //RK3 SSP TVD
    alpha1  = 1.00*dt;
    alpha21 = 0.25*dt;
    alpha22 = 0.25*dt;

    beta1 = 1.0/6.0*dt;
    beta2 = 1.0/6.0*dt;
    beta3 = 2.0/3.0*dt;
            
    
    //save Q
    for(i=0;i<*dof;++i){
        Q_stage[i] = Q[i];
    }
    
    
    /* stage 1 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k1);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k1);
    for(i=0;i<*dof;++i){
        Q[i] = Q_stage[i] + alpha1*k1[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    /* stage 2 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k2);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k2);
    for(i=0;i<*dof;++i){
        Q[i] = Q_stage[i] + alpha21*k1[i] + alpha22*k2[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    /* stage 3 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k3);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k3);
    for(i=0;i<*dof;++i){
        Q[i] = Q_stage[i] + beta1*k1[i] + beta2*k2[i] + beta3*k3[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    MPI_Allreduce(&resnorm_in,resnorm,1,MPI_DOUBLE,MPI_SUM,mpicomm_ext);
    *resnorm = sqrt(*resnorm);
    
}


void rk_4(  int* dim,
            int* unstructured,
            int* nfields,
            int* nelem_subgrid,
            int* ncell,
            int* nghost,
            int* nface,
            int* dof,
            int* cell_info,
            int* face_info,
            int* nsend,
            int* nrecv,
            int* send,
            int* recv,
            int* sizes,
            double* geom,
            double* Q,
            double* Q_stage,
            double* k1,
            double* k2,
            double* k3,
            double* k4,
            double dt,
            MPI_Comm mpicomm_ext,
            double* sendbuf,
            MPI_Request* request,
            int* counts,
            double* resnorm,
            double* compute_time,
            double* residual_time){
    
    int i;
    double t1,t2;
    double alpha1;
    double alpha2;
    double alpha3;
    double beta1;
    double resnorm_in;
    
    
    //RK2 SSP TVD
    alpha1 = 0.5*dt;
    alpha2 = 0.5*dt;
    alpha3 = 1.0*dt;
    beta1  = 1.0/6.0*dt;
            
    
    //save Q
    for(i=0;i<*dof;++i){
        Q_stage[i] = Q[i];
    }
    
    
    /* stage 1 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k1);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k1);
    for(i=0;i<*dof;++i){
        Q[i] = Q_stage[i] + alpha1*k1[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    /* stage 2 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k2);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k2);
    for(i=0;i<*dof;++i){
        Q[i] = Q_stage[i] + alpha2*k2[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    /* stage 3 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k3);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k3);
    for(i=0;i<*dof;++i){
        Q[i] = Q_stage[i] + alpha3*k3[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    /* stage 4 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k4);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k4);
    for(i=0;i<*dof;++i){
        Q[i] = Q_stage[i] + beta1*(k1[i]+ 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    MPI_Allreduce(&resnorm_in,resnorm,1,MPI_DOUBLE,MPI_SUM,mpicomm_ext);
    *resnorm = sqrt(*resnorm);
    
}


void lserk_45(  int* dim,
                int* unstructured,
                int* nfields,
                int* nelem_subgrid,
                int* ncell,
                int* nghost,
                int* nface,
                int* dof,
                int* cell_info,
                int* face_info,
                int* nsend,
                int* nrecv,
                int* send,
                int* recv,
                int* sizes,
                double* geom,
                double* Q,
                double* Q_stage,
                double* k1,
                double dt,
                MPI_Comm mpicomm_ext,
                double* sendbuf,
                MPI_Request* request,
                int* counts,
                double* resnorm,
                double* compute_time,
                double* residual_time){
    
    int i;
    double t1,t2;
    double       rk4a2,rk4a3,rk4a4,rk4a5;
    double rk4b1,rk4b2,rk4b3,rk4b4,rk4b5;
    double resnorm_in;
    
        
    //LSERK45
    rk4a2 = -567301805773.0  / 1357537059087.0;
    rk4a3 = -2404267990393.0 / 2016746695238.0;
    rk4a4 = -3550918686646.0 / 2091501179385.0;
    rk4a5 = -1275806237668.0 /  842570457699.0;

    rk4b1 = 1432997174477.0 /  9575080441755.0;
    rk4b2 = 5161836677717.0 / 13612068292357.0;
    rk4b3 = 1720146321549.0 /  2090206949498.0;
    rk4b4 = 3134564353537.0 /  4481467310338.0;
    rk4b5 = 2277821191437.0 / 14882151754819.0;
    
    
    
    /* stage 1 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k1);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k1);
    for(i=0;i<*dof;++i){
        Q_stage[i] = dt*k1[i];
        Q[i]       = Q[i] + rk4b1*Q_stage[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    /* stage 2 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k1);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k1);
    for(i=0;i<*dof;++i){
        Q_stage[i] = Q_stage[i]*rk4a2 + dt*k1[i];
        Q[i]       = Q[i] + rk4b2*Q_stage[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    /* stage 3 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k1);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k1);
    for(i=0;i<*dof;++i){
        Q_stage[i] = Q_stage[i]*rk4a3 + dt*k1[i];
        Q[i]       = Q[i] + rk4b3*Q_stage[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    /* stage 4 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k1);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k1);
    for(i=0;i<*dof;++i){
        Q_stage[i] = Q_stage[i]*rk4a4 + dt*k1[i];
        Q[i]       = Q[i] + rk4b4*Q_stage[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    /* stage 5 */
    physics_mpi_real_to_ghost(*nsend,*nrecv,*dof,send,recv,sizes,Q,mpicomm_ext,sendbuf,request,counts);
    
    t1 = MPI_Wtime();
    physics_rhs(dim,unstructured,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,Q,k1);
    t2 = MPI_Wtime();
    *residual_time += t2 - t1;
    
    resnorm_in = physics_norm_L2_residual(*dim,*nfields,*nelem_subgrid,*ncell,cell_info,k1);
    for(i=0;i<*dof;++i){
        Q_stage[i] = Q_stage[i]*rk4a5 + dt*k1[i];
        Q[i]       = Q[i] + rk4b5*Q_stage[i];
    }
    t2 = MPI_Wtime();
    *compute_time += t2 - t1;
    
    
    MPI_Allreduce(&resnorm_in,resnorm,1,MPI_DOUBLE,MPI_SUM,mpicomm_ext);
    *resnorm = sqrt(*resnorm);
    
}
