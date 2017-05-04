#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/* auxiliary routines, found at end of program */

char **allocchar2darray(int n, int m);
void freechar2darray(char **a);
void printarray(char **data, int n, int m);
void rowcol(int rank, const int blocks[2], int *row, int *col);
int isLastRow(int row, const int blocks[2]);
int isLastCol(int col, const int blocks[2]);
int typeIdx(int row, int col, const int blocks[2]);

/* first method - alltoallw */
void alltoall(const int myrow, const int mycol, const int rank, const int size, 
                    const int blocks[2], const int blocksize, const int globalsizes[2], const int localsizes[2],
                    const char *const globalptr,  char **localdata) {
    /*
 *      * get send and recieve counts ready for alltoallw call. 
 *           * everyone will be recieving just one block from proc 0; 
 *                * most procs will be sending nothing to anyone. 
 *                     */
    int sendcounts[ size ];
    int senddispls[ size ];
    MPI_Datatype sendtypes[size];
    int recvcounts[ size ];
    int recvdispls[ size ];
    MPI_Datatype recvtypes[size];

    for (int proc=0; proc<size; proc++) {
        recvcounts[proc] = 0;
        recvdispls[proc] = 0;
        recvtypes[proc] = MPI_CHAR;

        sendcounts[proc] = 0;
        senddispls[proc] = 0;
        sendtypes[proc] = MPI_CHAR;
    }
    recvcounts[0] = localsizes[0]*localsizes[1];
    recvdispls[0] = 0;


    /* The originating process needs to allocate and fill the source array,
 *      * and then define types defining the array chunks to send, and 
 *           * fill out senddispls, sendcounts (1) and sendtypes.
 *                */
    if (rank == 0) {
        /* 4 types of blocks - 
 *          * blocksize*blocksize, blocksize+1*blocksize, blocksize*blocksize+1, blocksize+1*blocksize+1
 *                   */
        MPI_Datatype blocktypes[4];
        int subsizes[2];
        int starts[2] = {0,0};
        for (int i=0; i<2; i++) {
           subsizes[0] = blocksize+i;
           for (int j=0; j<2; j++) {
               subsizes[1] = blocksize+j;
               MPI_Type_create_subarray(2, globalsizes, subsizes, starts, MPI_ORDER_C, MPI_CHAR, &blocktypes[2*i+j]);
               MPI_Type_commit(&blocktypes[2*i+j]);
           }
        }

        /* now figure out the displacement and type of each processor's data */
        for (int proc=0; proc<size; proc++) {
            int row, col;
            rowcol(proc, blocks, &row, &col);

            sendcounts[proc] = 1;
            senddispls[proc] = (row*blocksize*globalsizes[1] + col*blocksize)*sizeof(char);

            int idx = typeIdx(row, col, blocks);
            sendtypes[proc] = blocktypes[idx];
        }
    }

    MPI_Alltoallw(globalptr, sendcounts, senddispls, sendtypes,
                  &(localdata[0][0]), recvcounts, recvdispls, recvtypes, 
                  MPI_COMM_WORLD);
}


/* second  method: distribute almost all data using colums of size blocksize, 
 *  * then clean up the last row with another scatterv */

void twophasevecs(const int myrow, const int mycol, const int rank, const int size, 
                    const int blocks[2], const int blocksize, const int globalsizes[2], const int localsizes[2],
                    const char *const globalptr,  char **localdata) {
    int sendcounts[ size ];
    int senddispls[ size ];
    int recvcounts;

    for (int proc=0; proc<size; proc++) {
        sendcounts[proc] = 0;
        senddispls[proc] = 0;
    }

    /* We're going to be operating mostly in units of a single column of a "normal" sized block.
 *      * There will need to be two vectors describing these columns; one in the context of the
 *           * global array, and one in the local results.
 *                */
    MPI_Datatype vec, localvec;
    MPI_Type_vector(blocksize, 1, localsizes[1], MPI_CHAR, &localvec);
    MPI_Type_create_resized(localvec, 0, sizeof(char), &localvec);
    MPI_Type_commit(&localvec);

    MPI_Type_vector(blocksize, 1, globalsizes[1], MPI_CHAR, &vec);
    MPI_Type_create_resized(vec, 0, sizeof(char), &vec);
    MPI_Type_commit(&vec);

    /* The originating process needs to allocate and fill the source array,
 *      * and then define types defining the array chunks to send, and 
 *           * fill out senddispls, sendcounts (1) and sendtypes.
 *                */
    if (rank == 0) {
        /* create the vector type which will send one column of a "normal" sized-block */
        /* then all processors except those in the last row need to get blocksize*vec or (blocksize+1)*vec */
        /* will still have to do something to tidy up the last row of values */
        /* we need to make the type have extent of 1 char for scattering */
        for (int proc=0; proc<size; proc++) {
            int row, col;
            rowcol(proc, blocks, &row, &col);

            sendcounts[proc] = isLastCol(col, blocks) ? blocksize+1 : blocksize;
            senddispls[proc] = (row*blocksize*globalsizes[1] + col*blocksize);
        }
    }

    recvcounts = localsizes[1];
    MPI_Scatterv(globalptr, sendcounts, senddispls, vec,
                  &(localdata[0][0]), recvcounts, localvec, 0, MPI_COMM_WORLD);

    MPI_Type_free(&localvec);
    if (rank == 0)
        MPI_Type_free(&vec);

    /* now we need to do one more scatter, scattering just the last row of data 
 *      * just to the processors on the last row.
 *           * Here we recompute the sendcounts
 *                */
    if (rank == 0) {
        for (int proc=0; proc<size; proc++) {
            int row, col;
            rowcol(proc, blocks, &row, &col);
            sendcounts[proc] = 0;
            senddispls[proc] = 0;

            if ( isLastRow(row,blocks) ) {
                sendcounts[proc] = blocksize;
                senddispls[proc] = (globalsizes[0]-1)*globalsizes[1]+col*blocksize;
                if ( isLastCol(col,blocks) ) 
                    sendcounts[proc] += 1;
            }
        }
    }

    recvcounts = 0;
    if ( isLastRow(myrow, blocks) ) {
        recvcounts = blocksize;
        if ( isLastCol(mycol, blocks) )
            recvcounts++;
    }
    MPI_Scatterv(globalptr, sendcounts, senddispls, MPI_CHAR,
                  &(localdata[blocksize][0]), recvcounts, MPI_CHAR, 0, MPI_COMM_WORLD);
}
/* third method: first distribute rows, then columns, each with a single scatterv */

void twophaseRowCol(const int myrow, const int mycol, const int rank, const int size, 
                    const int blocks[2], const int blocksize, const int globalsizes[2], const int localsizes[2],
                    const char *const globalptr,  char **localdata) {
    char **rowdata ;

    /* create communicators which have processors with the same row or column in them*/
    MPI_Comm colComm, rowComm;
    MPI_Comm_split(MPI_COMM_WORLD, myrow, rank, &rowComm);
    MPI_Comm_split(MPI_COMM_WORLD, mycol, rank, &colComm);

    /* first, scatter the array by rows, with the processor in column 0 corresponding to each row
 *      * receiving the data */
    if (mycol == 0) {
        int sendcounts[ blocks[0] ];
        int senddispls[ blocks[0] ];
        senddispls[0] = 0;

        for (int row=0; row<blocks[0]; row++) {
            /* each processor gets blocksize rows, each of size globalsizes[1]... */
            sendcounts[row] = blocksize*globalsizes[1];
            if (row > 0) 
                senddispls[row] = senddispls[row-1] + sendcounts[row-1];
        }
        /* the last processor gets one more */
        sendcounts[blocks[0]-1] += globalsizes[1];

        /* allocate my rowdata */
        rowdata = allocchar2darray( sendcounts[myrow], globalsizes[1] );

        /* perform the scatter of rows */
        MPI_Scatterv(globalptr, sendcounts, senddispls, MPI_CHAR,
                      &(rowdata[0][0]), sendcounts[myrow], MPI_CHAR, 0, colComm);

    }

    /* Now, within each row of processors, we can scatter the columns.  
 *      * We can do this as we did in the previous example; create a vector
 *           * (and localvector) type and scatter accordingly */
    int locnrows = blocksize;
    if ( isLastRow(myrow, blocks) )
        locnrows++;

    MPI_Datatype vec, localvec;
    MPI_Type_vector(locnrows, 1, globalsizes[1], MPI_CHAR, &vec);
    MPI_Type_create_resized(vec, 0, sizeof(char), &vec);
    MPI_Type_commit(&vec);

    MPI_Type_vector(locnrows, 1, localsizes[1], MPI_CHAR, &localvec);
    MPI_Type_create_resized(localvec, 0, sizeof(char), &localvec);
    MPI_Type_commit(&localvec);

    int sendcounts[ blocks[1] ];
    int senddispls[ blocks[1] ];
    if (mycol == 0) {
        for (int col=0; col<blocks[1]; col++) {
            sendcounts[col] = isLastCol(col, blocks) ? blocksize+1 : blocksize;
            senddispls[col] = col*blocksize;
        }
    }
    char *rowptr = (mycol == 0) ? &(rowdata[0][0]) : NULL;

    MPI_Scatterv(rowptr, sendcounts, senddispls, vec,
                  &(localdata[0][0]), sendcounts[mycol], localvec, 0, rowComm);

    MPI_Type_free(&localvec);
    MPI_Type_free(&vec);

    if (mycol == 0) 
        freechar2darray(rowdata);

    MPI_Comm_free(&rowComm);
    MPI_Comm_free(&colComm);
}

int main(int argc, char **argv) {

    int rank, size;
    int blocks[2] = {0,0};
    const int blocksize=3;
    int globalsizes[2], localsizes[2];
    char **globaldata;
    char *globalptr = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0 && argc < 2) {
        fprintf(stderr,"Usage: %s method\n   Where method is one of: alltoall, twophasevecs, twophaserowcol\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    /* calculate sizes for a 2d grid of processors */
    MPI_Dims_create(size, 2, blocks);

    int myrow, mycol;
    rowcol(rank, blocks, &myrow, &mycol);

    /* create array sizes so that last block has 1 too many rows/cols */
    globalsizes[0] = blocks[0]*blocksize+1;  
    globalsizes[1] = blocks[1]*blocksize+1;
    if (rank == 0) {
        globaldata = allocchar2darray(globalsizes[0], globalsizes[1]);
        globalptr = &(globaldata[0][0]);
        for (int i=0; i<globalsizes[0]; i++) 
            for (int j=0; j<globalsizes[1]; j++)
                globaldata[i][j] = 'a'+(i*globalsizes[1] + j)%26;

        printf("Global array: \n");
        printarray(globaldata, globalsizes[0], globalsizes[1]);
    }

    /* the local chunk we'll be receiving */
    localsizes[0] = blocksize; localsizes[1] = blocksize;
    if ( isLastRow(myrow,blocks)) localsizes[0]++;
    if ( isLastCol(mycol,blocks)) localsizes[1]++;
    char **localdata = allocchar2darray(localsizes[0],localsizes[1]);

    if (!strcasecmp(argv[1], "alltoall")) {
        if (rank == 0) printf("Method - alltoall\n");
        alltoall(myrow, mycol, rank, size, blocks, blocksize, globalsizes, localsizes, globalptr, localdata);
    } else if (!strcasecmp(argv[1],"twophasevecs")) {
        if (rank == 0) printf("Method - two phase, vectors, then cleanup\n");
        twophasevecs(myrow, mycol, rank, size, blocks, blocksize, globalsizes, localsizes, globalptr, localdata);
    } else {
        if (rank == 0) printf("Method - two phase - row, cols\n");
        twophaseRowCol(myrow, mycol, rank, size, blocks, blocksize, globalsizes, localsizes, globalptr, localdata);
    }

    for (int proc=0; proc<size; proc++) {
        if (proc == rank) {
            printf("\nRank %d:\n", proc);
            printarray(localdata, localsizes[0], localsizes[1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);            
    }

    freechar2darray(localdata);
    if (rank == 0) 
        freechar2darray(globaldata);

    MPI_Finalize();

    return 0;
}

char **allocchar2darray(int n, int m) {
    char **ptrs = malloc(n*sizeof(char *));
    ptrs[0] = malloc(n*m*sizeof(char));
    for (int i=0; i<n*m; i++)
        ptrs[0][i]='.';

    for (int i=1; i<n; i++) 
        ptrs[i] = ptrs[i-1] + m;

    return ptrs;
}

void freechar2darray(char **a) {
    free(a[0]);
    free(a);
}

void printarray(char **data, int n, int m) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) 
            putchar(data[i][j]);
        putchar('\n');
    }
}

void rowcol(int rank, const int blocks[2], int *row, int *col) {
    *row = rank/blocks[1];
    *col = rank % blocks[1];
}

int isLastRow(int row, const int blocks[2]) {
    return (row == blocks[0]-1);
}

int isLastCol(int col, const int blocks[2]) {
    return (col == blocks[1]-1);
}

int typeIdx(int row, int col, const int blocks[2]) {
    int lastrow = (row == blocks[0]-1);
    int lastcol = (col == blocks[1]-1);

    return lastrow*2 + lastcol;
}
