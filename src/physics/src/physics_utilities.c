/**
 * \file   physics_utilities.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Utility functions for the PHYSICS code module.
 * 
 * Created on August 21, 2018, 2:10 PM
 */

/* header files */
#include "physics_utilities.h"


void utilities_write_physics_message(int mpi_rank){
    
    if(mpi_rank==0){
     
        printf("\n");
        printf("+==========================================+\n");
        printf("|              Physics Kernel              |\n");
        printf("+==========================================+\n");
        
    }

}


int utilities_find_keyword_integer(char* filename, char *keyword, int *integer){
    FILE *fp;
    char buff[BUFFER_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
            *integer = atoi(&buff[length]);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_two_integers(char* filename, char *keyword, 
                                        int *integer1, 
                                        int *integer2){
    FILE *fp;
    char buff[BUFFER_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
      
          sscanf(&buff[length],"%d %d", integer1,integer2);
          sscanf(&buff[length],"%d, %d",integer1,integer2);
          sscanf(&buff[length],"%d,%d",integer1,integer2);
          
          fclose(fp);
          return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_three_integers(char* filename, char *keyword, 
                                            int *integer1, 
                                            int *integer2, 
                                            int *integer3){
    FILE *fp;
    char buff[BUFFER_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
      
            sscanf(&buff[length],"%d %d %d",  integer1,integer2,integer3);
            sscanf(&buff[length],"%d, %d, %d",integer1,integer2,integer3);
            sscanf(&buff[length],"%d, %d %d", integer1,integer2,integer3);
            sscanf(&buff[length],"%d %d, %d", integer1,integer2,integer3);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_four_integers(char* filename, char *keyword, 
                                         int *integer1, 
                                         int *integer2, 
                                         int *integer3,
                                         int *integer4){
    FILE *fp;
    char buff[BUFFER_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
      
            sscanf(&buff[length],"%d %d %d %d",  integer1,integer2,integer3,integer4);
            sscanf(&buff[length],"%d, %d, %d, %d",integer1,integer2,integer3,integer4);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_five_integers(char* filename, char *keyword, 
                                         int *integer1, 
                                         int *integer2, 
                                         int *integer3,
                                         int *integer4,
                                         int *integer5){
    FILE *fp;
    char buff[BUFFER_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
      
            sscanf(&buff[length],"%d %d %d %d %d",  integer1,integer2,integer3,integer4,integer5);
            sscanf(&buff[length],"%d, %d, %d, %d, %d",integer1,integer2,integer3,integer4,integer5);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_double(char* filename, char *keyword, double *dbl){
    FILE *fp;
    char buff[BUFFER_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
            *dbl = atof(&buff[length]);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_two_doubles(char* filename, char *keyword, 
                                       double *dbl1, 
                                       double *dbl2){
    FILE *fp;
    char buff[BUFFER_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
            sscanf(&buff[length],"%le %le",dbl1,dbl2);
            sscanf(&buff[length],"%le, %le",dbl1,dbl2);
            sscanf(&buff[length],"%le,%le",dbl1,dbl2);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_three_doubles(char* filename, char *keyword, 
                                            double *dbl1, 
                                            double *dbl2, 
                                            double *dbl3){
    FILE *fp;
    char buff[BUFFER_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
            sscanf(&buff[length],"%le %le %le",  dbl1,dbl2,dbl3);
            sscanf(&buff[length],"%le, %le, %le",dbl1,dbl2,dbl3);
            sscanf(&buff[length],"%le %le, %le", dbl1,dbl2,dbl3);
            sscanf(&buff[length],"%le, %le %le", dbl1,dbl2,dbl3);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_four_doubles(char* filename, char *keyword, 
                                        double *dbl1, 
                                        double *dbl2, 
                                        double *dbl3,
                                        double *dbl4){
    FILE *fp;
    char buff[BUFFER_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
            sscanf(&buff[length],"%le %le %le %le",   dbl1,dbl2,dbl3,dbl4);
            sscanf(&buff[length],"%le, %le, %le, %le",dbl1,dbl2,dbl3,dbl4);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_five_doubles(char* filename, char *keyword, 
                                        double *dbl1, 
                                        double *dbl2, 
                                        double *dbl3,
                                        double *dbl4,
                                        double *dbl5){
    FILE *fp;
    char buff[BUFFER_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
            sscanf(&buff[length],"%le %le %le %le %le",    dbl1,dbl2,dbl3,dbl4,dbl5);
            sscanf(&buff[length],"%le, %le, %le, %le, %le",dbl1,dbl2,dbl3,dbl4,dbl5);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_string(char *filename, char *keyword, 
                                    char *string){
    FILE *fp;
    char buff[BUFFER_SIZE];
    unsigned long i,length;

    fp = fopen(filename,"r");
    length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
      
            int ind=0;
            for(i=length;i<BUFFER_SIZE;i++){
                if(buff[i] != ' '){
                    string[ind] = buff[i];
                    ind = ind+1;
                }
            }

            string[strcspn(string, "\n")] = 0;

            if(strlen(string) >= BUFFER_SIZE){
                printf("[ AMR ] string length is greater than allowable");
                exit(0);
            }

            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}
