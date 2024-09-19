/**
 * \file   physics_runge_kutta_interface.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Explicit Runge-Kutta time stepping interface for solving
 *        time-accurate physics using Method-of-Lines.
 *
 * Explicit Runge-Kutta methods available:
 *      * ERK1: one stage, first-order accuracy
 *      * ERK2: two stage, second-order accuracy
 *      * ERK3: three stage, third-order accuracy
 *      * ERK4: four stage, fourth-order accuracy
 *      * LSERK45: five stage, fourth-order accuracy, low storage
 *
 * Created on August 24, 2018, 2:53 PM
 */

/* header files */
#include "physics_runge_kutta_interface.h"


void physics_runge_kutta_interface( int* dim,
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
                                    int* istep,
                                    int* regrid,
                                    int* visualize,
                                    int* checkpoint,
                                    MPI_Comm mpicomm_ext,
                                    int* counts,
                                    double* sendbuf,
                                    MPI_Request* request,
                                    int* step_count,
                                    double* t1,
                                    double* t2,
                                    double* compute_time,
                                    double* residual_time,
                                    double* residual_norm){


    double *k1;
    double *k2;
    double *k3;
    double *k4;
    double *Q_stage;


    double dt;
    double dt_local;
    double time1;
    double time2;

    int go_solver;


    *step_count = 0;
    *compute_time = 0.0;
    *residual_time = 0.0;
    *residual_norm = -1.0;



    go_solver = 1;
    switch(d_physics_simulation.time_scheme){
        case(1): //RK1

            /* allocate memory */
            k1 = (double *) malloc(sizeof(double)*(*dof));

            *t1 = MPI_Wtime();
            while(go_solver){

                time1 = MPI_Wtime();
                    dt = d_physics_simulation.dt;
                    physics_time_step(dim,nfields,nelem_subgrid,
                                      ncell,cell_info,geom,Q,&dt);
                time2 = MPI_Wtime();
                *compute_time += (time2 - time1);

                dt_local = dt;
                MPI_Allreduce(&dt_local,&dt,1,MPI_DOUBLE,MPI_MIN,mpicomm_ext);

                rk_1(   dim,unstructured,nfields,nelem_subgrid,ncell,nghost,nface,
                        dof,cell_info,face_info,nsend,nrecv,send,recv,sizes,geom,
                        Q,k1,dt,mpicomm_ext,
                        sendbuf,request,counts,
                        residual_norm,compute_time,residual_time);


                *istep += 1;
                *step_count += 1;
                if(*istep>=d_physics_simulation.time_steps){
                    go_solver = 0;
                }
                if(*istep%d_physics_simulation.regrid_interval==0){
                    go_solver = 0;
                    *regrid = 1;
                }
                if(*istep%d_physics_simulation.checkpoint_interval==0){
                    go_solver = 0;
                    *checkpoint = 1;
                }
                if(*istep%d_physics_simulation.visualization_interval==0){
                    go_solver = 0;
                    *visualize = 1;
                }
            }
            *t2 = MPI_Wtime();


            /* free memory */
            free(k1);
            break;

        case(2): //RK2

            /* allocate memory */
            k1 = (double *) malloc(sizeof(double)*(*dof));
            k2 = (double *) malloc(sizeof(double)*(*dof));
            Q_stage = (double *) malloc(sizeof(double)*(*dof));


            *t1 = MPI_Wtime();
            while(go_solver){

                time1 = MPI_Wtime();
                    dt = d_physics_simulation.dt;
                //    physics_time_step(dim,nfields,nelem_subgrid,
                //                      ncell,cell_info,geom,Q,&dt);
                time2 = MPI_Wtime();
                *compute_time += (time2 - time1);

                dt_local = dt;
                MPI_Allreduce(&dt_local,&dt,1,MPI_DOUBLE,MPI_MIN,mpicomm_ext);
		//printf("Time step: %e\n",dt);

                rk_2(   dim,unstructured,nfields,nelem_subgrid,ncell,nghost,nface,
                        dof,cell_info,face_info,nsend,nrecv,send,recv,sizes,geom,
                        Q,Q_stage,k1,k2,dt,mpicomm_ext,
                        sendbuf,request,counts,
                        residual_norm,compute_time,residual_time);

                *istep += 1;
                *step_count += 1;
                if(*istep>=d_physics_simulation.time_steps){
                    go_solver = 0;
                }
                if(*istep%d_physics_simulation.regrid_interval==0){
                    go_solver = 0;
                    *regrid = 1;
                }
                if(*istep%d_physics_simulation.checkpoint_interval==0){
                    go_solver = 0;
                    *checkpoint = 1;
                }
                if(*istep%d_physics_simulation.visualization_interval==0){
                    go_solver = 0;
                    *visualize = 1;
                }
            }
            *t2 = MPI_Wtime();


            /* free memory */
            free(k1);
            free(k2);
            free(Q_stage);
            break;

        case(3): //RK3

            /* allocate memory */
            k1 = (double *) malloc(sizeof(double)*(*dof));
            k2 = (double *) malloc(sizeof(double)*(*dof));
            k3 = (double *) malloc(sizeof(double)*(*dof));
            Q_stage = (double *) malloc(sizeof(double)*(*dof));


            *t1 = MPI_Wtime();
            while(go_solver){

                time1 = MPI_Wtime();
                    dt = d_physics_simulation.dt;
                    physics_time_step(dim,nfields,nelem_subgrid,
                                      ncell,cell_info,geom,Q,&dt);
                time2 = MPI_Wtime();
                *compute_time += (time2 - time1);

                dt_local = dt;
                MPI_Allreduce(&dt_local,&dt,1,MPI_DOUBLE,MPI_MIN,mpicomm_ext);

                rk_3(   dim,unstructured,nfields,nelem_subgrid,ncell,nghost,nface,
                        dof,cell_info,face_info,nsend,nrecv,send,recv,sizes,geom,
                        Q,Q_stage,k1,k2,k3,dt,mpicomm_ext,
                        sendbuf,request,counts,
                        residual_norm,compute_time,residual_time);

                *istep += 1;
                *step_count += 1;
                if(*istep>=d_physics_simulation.time_steps){
                    go_solver = 0;
                }
                if(*istep%d_physics_simulation.regrid_interval==0){
                    go_solver = 0;
                    *regrid = 1;
                }
                if(*istep%d_physics_simulation.checkpoint_interval==0){
                    go_solver = 0;
                    *checkpoint = 1;
                }
                if(*istep%d_physics_simulation.visualization_interval==0){
                    go_solver = 0;
                    *visualize = 1;
                }
            }
            *t2 = MPI_Wtime();


            /* free memory */
            free(k1);
            free(k2);
            free(k3);
            free(Q_stage);
            break;

        case(4): //RK4

            /* allocate memory */
            k1 = (double *) malloc(sizeof(double)*(*dof));
            k2 = (double *) malloc(sizeof(double)*(*dof));
            k3 = (double *) malloc(sizeof(double)*(*dof));
            k4 = (double *) malloc(sizeof(double)*(*dof));
            Q_stage = (double *) malloc(sizeof(double)*(*dof));


            *t1 = MPI_Wtime();
            while(go_solver){

                time1 = MPI_Wtime();
                    dt = d_physics_simulation.dt;
                    physics_time_step(dim,nfields,nelem_subgrid,
                                      ncell,cell_info,geom,Q,&dt);
                time2 = MPI_Wtime();
                *compute_time += (time2 - time1);

                dt_local = dt;
                MPI_Allreduce(&dt_local,&dt,1,MPI_DOUBLE,MPI_MIN,mpicomm_ext);

                rk_4(   dim,unstructured,nfields,nelem_subgrid,ncell,nghost,nface,
                        dof,cell_info,face_info,nsend,nrecv,send,recv,sizes,geom,
                        Q,Q_stage,k1,k2,k3,k4,dt,mpicomm_ext,
                        sendbuf,request,counts,
                        residual_norm,compute_time,residual_time);

                *istep += 1;
                *step_count += 1;
                if(*istep>=d_physics_simulation.time_steps){
                    go_solver = 0;
                }
                if(*istep%d_physics_simulation.regrid_interval==0){
                    go_solver = 0;
                    *regrid = 1;
                }
                if(*istep%d_physics_simulation.checkpoint_interval==0){
                    go_solver = 0;
                    *checkpoint = 1;
                }
                if(*istep%d_physics_simulation.visualization_interval==0){
                    go_solver = 0;
                    *visualize = 1;
                }
            }
            *t2 = MPI_Wtime();


            /* free memory */
            free(k1);
            free(k2);
            free(k3);
            free(k4);
            free(Q_stage);
            break;

        case(5): //LSERK45

            /* allocate memory */
            k1 = (double *) malloc(sizeof(double)*(*dof));
            Q_stage = (double *) malloc(sizeof(double)*(*dof));

            *t1 = MPI_Wtime();
            while(go_solver){

                time1 = MPI_Wtime();
                    dt = d_physics_simulation.dt;
                    physics_time_step(dim,nfields,nelem_subgrid,
                                      ncell,cell_info,geom,Q,&dt);
                time2 = MPI_Wtime();
                *compute_time += (time2 - time1);

                dt_local = dt;
                MPI_Allreduce(&dt_local,&dt,1,MPI_DOUBLE,MPI_MIN,mpicomm_ext);

                lserk_45(   dim,unstructured,nfields,nelem_subgrid,ncell,nghost,nface,
                            dof,cell_info,face_info,nsend,nrecv,send,recv,sizes,geom,
                            Q,Q_stage,k1,dt,mpicomm_ext,
                            sendbuf,request,counts,
                            residual_norm,compute_time,residual_time);

                *istep += 1;
                *step_count += 1;
                if(*istep>=d_physics_simulation.time_steps){
                    go_solver = 0;
                }
                if(*istep%d_physics_simulation.regrid_interval==0){
                    go_solver = 0;
                    *regrid = 1;
                }
                if(*istep%d_physics_simulation.checkpoint_interval==0){
                    go_solver = 0;
                    *checkpoint = 1;
                }
                if(*istep%d_physics_simulation.visualization_interval==0){
                    go_solver = 0;
                    *visualize = 1;
                }
            }
            *t2 = MPI_Wtime();


            /* free memory */
            free(k1);
            free(Q_stage);
            break;

    }

    if(isnan(*residual_norm)) {
        printf("[solve] nan detected in residual aborting\n");
        MPI_Abort(mpicomm_ext,3);
    }

}

