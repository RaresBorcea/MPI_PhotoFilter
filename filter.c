#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<stddef.h>

#define MASTER      0       // Master proc
#define NAMELEN     30      // Max file name length
#define FILTNO      20      // Max no of filters to apply
#define LINELENGHT  50      // Max chars per line for input
#define VALLEN      10      // Max chars for image params
#define TYPELEN     3       // Image type length
#define KERSIZE     3       // Filters kernel size
#define MAXPROC     48      // Max no of procs
#define NOOFFILT    5       // No of possible filters
#define BWIMG       "P5"    // B&W image
#define COLORIMG    "P6"    // Color image

// RGB structure for PNM PPM file/Color image
typedef struct matrix {
    unsigned char red; 
    unsigned char green; 
    unsigned char blue;
} __attribute__((packed)) RGBMatrix;

// Get cli arguments and store them accordingly
void getArgs(int argc, char **argv, char *file_in, char *file_out, int *filtToApply, int *filtSize) {
    if (argc < 4) {
        printf("Not enough paramters: mpirun -np N ./tema3 image_in.pnm image_out.pnm \
            filter1 filter2 ... filterX\n");
        exit(1);
    }

    // Store input file names
    strcpy(file_in, argv[1]);
    strcpy(file_out, argv[2]);

    // Store fitlers to be applied
    // Number-code each filter type:
    // 1=smooth, 2=blur, 3=sharpen, 4=mean, 5=emboss
    int j = 0;
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "smooth") == 0) {
            filtToApply[j++] = 1;
        } if (strcmp(argv[i], "blur") == 0) {
            filtToApply[j++] = 2;
        } else if (strcmp(argv[i], "sharpen") == 0) {
            filtToApply[j++] = 3;
        } else if (strcmp(argv[i], "mean") == 0) {
            filtToApply[j++] = 4;
        } else if (strcmp(argv[i], "emboss") == 0) {
            filtToApply[j++] = 5;
        }
    }

    *filtSize = j;
}

// Read input file
void readFile(char *file, char *type, int *width, int *height, int *vmax, 
    RGBMatrix ***rgbMatrix, unsigned char ***bwMatrix, int rank) {
    FILE *input;

    // Only MASTER proc reads the file
    if (rank == MASTER) {
        char line[LINELENGHT], widthchar[VALLEN], heightchar[VALLEN], vmaxchar[VALLEN];
        int i, j, k;

        input = fopen(file, "rb");
        if (!input) {
            printf("Couldn't open input file. Exiting.");
            exit(1);
        }

        // Read image type
        fread(type, sizeof(char), TYPELEN, input); 
        type[TYPELEN - 1] = 0;

        // Read comment line and discard it
        i = 0;
        do {
            fread(&line[i], sizeof(char), 1, input);
            i++;
        } while (line[i - 1] != '\n');

        // Read image width and height
        i = 0;
        do {
            fread(&line[i], sizeof(char), 1, input);
            i++;
        } while (line[i - 1] != '\n');
        i--;
        for (j = 0; isspace(line[j]) == 0; j++) {
            widthchar[j] = line[j];
        }
        widthchar[j] = 0;

        *width = atoi(widthchar);

        while (!isdigit(line[j])) {
            j++;
        }
        k = 0;
        while (j < i) {
            heightchar[k] = line[j];
            k++; j++;
        }
        heightchar[k] = 0;

        *height = atoi(heightchar);

        // Read max value
        i = 0;
        do {
            fread(&line[i], sizeof(char), 1, input);
            i++;
        } while (line[i - 1] != '\n');
        i--;
        j = 0;
        while (j < i) {
            vmaxchar[j] = line[j];
            j++;
        }
        vmaxchar[j] = 0;

        *vmax = atoi(vmaxchar);
    }
    
    // Broadcast image type, width and height
    MPI_Bcast(type, TYPELEN, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(height, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(width, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (strcmp(type, BWIMG) == 0) {
        // Allocate matrix for B&W image
        *bwMatrix = (unsigned char**) malloc ((*height + 2) * sizeof(unsigned char*));
        (*bwMatrix)[0] = (unsigned char*) calloc ((*height + 2) * (*width + 2), sizeof(unsigned char));
        for (int i = 1; i < *height + 2; i++) {
            (*bwMatrix)[i] = (*bwMatrix)[0] + i * (*width + 2);
        }
    } else if (strcmp(type, COLORIMG) == 0) {
        // Allocate matrix for Color image
        *rgbMatrix = (RGBMatrix**) malloc ((*height + 2) * sizeof(RGBMatrix*));
        (*rgbMatrix)[0] = (RGBMatrix*) calloc ((*height + 2) * (*width + 2), sizeof(RGBMatrix));
        for (int i = 1; i < *height + 2; i++) {
            (*rgbMatrix)[i] = (*rgbMatrix)[0] + i * (*width + 2);
        }
    }

    if (rank == MASTER) {
        // Read image content into matrix
        if (strcmp(type, BWIMG) == 0) {
            // Respect ZERO border, save starting from index 1
            for (int i = 1; i < *height + 1; i++) {
                fread((*bwMatrix)[i] + 1, sizeof(unsigned char), *width, input);
            }
        } else if (strcmp(type, COLORIMG) == 0) {
            for (int i = 1; i < *height + 1; i++) {
                fread((*rgbMatrix)[i] + 1, sizeof(RGBMatrix), *width, input);
            }
        }

        fclose(input);
    }
}

// Write output file
void saveFile(char *fileout, int width, int height, RGBMatrix **rgbMatrix, 
    unsigned char **bwMatrix, int vmax, char *type) {

    FILE *output = fopen(fileout, "wb");
    if (!output) {
        printf("Couldn't open output file. Exiting.");
        exit(1);
    } else {
        fprintf(output, "%s\n%d %d\n%d\n", type, width, height, vmax);

        if (strcmp(type, BWIMG) == 0) {
            for (int i = 1; i < height + 1; i++) {
                fwrite(bwMatrix[i] + 1, sizeof(unsigned char), width, output);
            }
        } else if (strcmp(type, COLORIMG) == 0) {
            for (int i = 1; i < height + 1; i++) {
                fwrite(rgbMatrix[i] + 1, sizeof(RGBMatrix), width, output);
            }
        }
        
        fclose(output);
    }
}

// Rotate matrix by 180 degree 
void rotateMatrix180(int size, float a[][size], float rotatedA[][size]) { 
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++){
            rotatedA[i][j] = a[size - i - 1][size - j - 1];
        }
    } 
}

// Apply given filter on given input matrix and return the result in given output matrix
void applyFilter(int startPos, int width, float rotatedFilter[][KERSIZE], unsigned char** bwMatrixIn, 
    unsigned char** bwMatrixOut, RGBMatrix** rgbMatrixIn, RGBMatrix** rgbMatrixOut, int rowsPerProcLoc) { 
    // Check input matrix type to determine sum calculation approach
    if (rgbMatrixIn == NULL) {
        for (int i = startPos; i < startPos + rowsPerProcLoc; i++) {
            for (int j = 1; j < width + 1; j++) {
                // For B&W image, use neighbour pixels values directly
                float sum = 0;
                int r, s, m, n;
                for (r = i - 1, m = 0; r <= i + 1 && m < KERSIZE; r++, m++) {
                    for (s = j - 1, n = 0; s <= j + 1 && n < KERSIZE; s++, n++) {
                        sum += rotatedFilter[m][n] * bwMatrixIn[r][s];
                    }
                }

                // Clamp the value in 0-255 space
                if (sum > 255) {
                    sum = 255;
                } else if (sum < 0) {
                    sum = 0;
                }

                // Cast to image matrix type
                unsigned char newValue = (unsigned char)sum;

                bwMatrixOut[i][j] = newValue;
            }
        }
    } else {
        for (int i = startPos; i < startPos + rowsPerProcLoc; i++) {
            for (int j = 1; j < width + 1; j++) {
                // For Color image, use neighbour pixels values per RGB channels
                float sumR = 0, sumG = 0, sumB = 0;
                int r, s, m, n;
                for (r = i - 1, m = 0; r <= i + 1 && m < KERSIZE; r++, m++) {
                    for (s = j - 1, n = 0; s <= j + 1 && n < KERSIZE; s++, n++) {
                        sumR += rotatedFilter[m][n] * rgbMatrixIn[r][s].red;
                        sumG += rotatedFilter[m][n] * rgbMatrixIn[r][s].green;
                        sumB += rotatedFilter[m][n] * rgbMatrixIn[r][s].blue;
                    }
                }

                // Clamp values in 0-255 space
                if (sumR > 255) {
                    sumR = 255;
                } else if (sumR < 0) {
                    sumR = 0;
                }

                if (sumG > 255) {
                    sumG = 255;
                } else if (sumG < 0) {
                    sumG = 0;
                }

                if (sumB > 255) {
                    sumB = 255;
                } else if (sumB < 0) {
                    sumB = 0;
                }

                // Cast to image matrix type
                unsigned char newValueR = (unsigned char)sumR;
                unsigned char newValueG = (unsigned char)sumG;
                unsigned char newValueB = (unsigned char)sumB;

                RGBMatrix newValue;
                newValue.red = newValueR;
                newValue.green = newValueG;
                newValue.blue = newValueB;

                rgbMatrixOut[i][j] = newValue;
            }
        }
    }
}

// Driver function
int main(int argc, char * argv[]) {
    char file_in[NAMELEN];
    char file_out[NAMELEN];
    int filtToApply[FILTNO];    // filters to be applied, in order
    char type[TYPELEN];         // file type
    int width, height, vmax;

    // Matrices used to store image
    unsigned char **bwMatrix1;  // primary matrix - B&W image
    unsigned char **bwMatrix2;  // secondary matrix - B&W image
    RGBMatrix **rgbMatrix1;     // primary matrix - Color image
    RGBMatrix **rgbMatrix2;     // secondary matrix - Color image

    // Used to calculate image rows distribution per proc
    int offsetLoc, rowsPerProcLoc, averageRows, restRows;
    int rowsPerProc[MAXPROC];
    int offset[MAXPROC];

    // Unrotated filters kernel matrices
    float smooth[KERSIZE][KERSIZE] = {  {(float)1/9, (float)1/9, (float)1/9},  
                                        {(float)1/9, (float)1/9, (float)1/9},  
                                        {(float)1/9, (float)1/9, (float)1/9} };

    float blur[KERSIZE][KERSIZE] = {    {(float)1/16, (float)2/16, (float)1/16},
                                        {(float)2/16, (float)4/16, (float)2/16},  
                                        {(float)1/16, (float)2/16, (float)1/16}  };

    float sharpen[KERSIZE][KERSIZE] = { {0, -(float)2/3, 0},
                                        {-(float)2/3, (float)11/3, -(float)2/3},  
                                        {0, -(float)2/3, 0}    };

    float mean[KERSIZE][KERSIZE] = {    {-1, -1, -1},
                                        {-1, 9, -1},  
                                        {-1, -1, -1}    };

    float emboss[KERSIZE][KERSIZE] = {  {0, 1, 0},
                                        {0, 0, 0},  
                                        {0, -1, 0}  };

    // Rotate filters kernel matrices
    float rotatedSmooth[KERSIZE][KERSIZE];
    float rotatedBlur[KERSIZE][KERSIZE];
    float rotatedSharpen[KERSIZE][KERSIZE];
    float rotatedMean[KERSIZE][KERSIZE];
    float rotatedEmboss[KERSIZE][KERSIZE];

    // Used to check rotation status, for efficiency
    int alreadyRotated[NOOFFILT + 1];
    int alreadyBcasted[NOOFFILT + 1];

    memset(alreadyRotated, 0, sizeof(alreadyRotated));
    memset(alreadyBcasted, 0, sizeof(alreadyBcasted));
    memset(filtToApply, 0, sizeof(filtToApply));

    // Used to interchange input<->result matrices after each step
    int isT2forResult = 1;

    int rank;
    int nProcesses;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

    // Create derived MPI data type for RGBMatrix
    const int nitems = 3;
    int blocklens[3] = {1, 1, 1};
    MPI_Datatype old_types[3] = {MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR};
    MPI_Datatype mpi_rgbmatrix;
    MPI_Aint offsets[3];

    offsets[0] = offsetof(RGBMatrix, red);
    offsets[1] = offsetof(RGBMatrix, green);
    offsets[2] = offsetof(RGBMatrix, blue);

    MPI_Type_create_struct(nitems, blocklens, offsets, old_types, &mpi_rgbmatrix);
    MPI_Type_commit(&mpi_rgbmatrix);

    // Get cli args and bcast filters to be applied
    int filtSize;
    if (rank == MASTER) { 
        getArgs(argc, argv, file_in, file_out, filtToApply, &filtSize);
    }
    MPI_Bcast(&filtSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filtToApply, filtSize, MPI_INT, 0, MPI_COMM_WORLD);

    readFile(file_in, type, &width, &height, &vmax, &rgbMatrix1, &bwMatrix1, rank);

    // Allocate secondary matrix according to image type
    if (strcmp(type, BWIMG) == 0) {
        bwMatrix2 = (unsigned char**) malloc ((height + 2) * sizeof(unsigned char*));
        bwMatrix2[0] = (unsigned char*) calloc ((height + 2) * (width + 2), sizeof(unsigned char));
        for (int i = 1; i < height + 2; i++) {
            bwMatrix2[i] = bwMatrix2[0] + i * (width + 2);
        }
    } else if (strcmp(type, COLORIMG) == 0) {
        rgbMatrix2 = (RGBMatrix**) malloc ((height + 2) * sizeof(RGBMatrix*));
        rgbMatrix2[0] = (RGBMatrix*) calloc ((height + 2) * (width + 2), sizeof(RGBMatrix));
        for (int i = 1; i < height + 2; i++) {
            rgbMatrix2[i] = rgbMatrix2[0] + i * (width + 2);
        }
    }

    // Rotate filters kernel matrices in MASTER proc
    if (rank == MASTER) {
        int i = 0;
        while (filtToApply[i] != 0) {
            switch (filtToApply[i]) { 
                case 1: if (alreadyRotated[1] == 0) {
                            rotateMatrix180(KERSIZE, smooth, rotatedSmooth);
                            // Rotate once, after first appearance in toApply
                            alreadyRotated[1] = 1;
                        }
                        break; 
                case 2: if (alreadyRotated[2] == 0) {
                            rotateMatrix180(KERSIZE, blur, rotatedBlur);
                            alreadyRotated[2] = 1;
                        }
                        break; 
                case 3: if (alreadyRotated[3] == 0) {
                            rotateMatrix180(KERSIZE, sharpen, rotatedSharpen);
                            alreadyRotated[3] = 1;
                        }
                        break; 
                case 4: if (alreadyRotated[4] == 0) {
                            rotateMatrix180(KERSIZE, mean, rotatedMean);
                            alreadyRotated[4] = 1;
                        }
                        break; 
                case 5: if (alreadyRotated[5] == 0) {
                            rotateMatrix180(KERSIZE, emboss, rotatedEmboss);
                            alreadyRotated[5] = 1;
                        }
                        break; 
            } 
            i++;
        }
    }

    // Bcast rotated filters kernel matrices
    int i = 0;
    while (filtToApply[i] != 0) {
        switch (filtToApply[i]) { 
            case 1: if (alreadyBcasted[1] == 0) {
                        MPI_Bcast(rotatedSmooth, KERSIZE * KERSIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);
                        // Bcast once, after first appearance in toApply
                        alreadyBcasted[1] = 1;
                    }
                    break; 
            case 2: if (alreadyBcasted[2] == 0) {
                        MPI_Bcast(rotatedBlur, KERSIZE * KERSIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);
                        alreadyBcasted[2] = 1;
                    }
                    break; 
            case 3: if (alreadyBcasted[3] == 0) {
                        MPI_Bcast(rotatedSharpen, KERSIZE * KERSIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);
                        alreadyBcasted[3] = 1;
                    }
                    break; 
            case 4: if (alreadyBcasted[4] == 0) {
                        MPI_Bcast(rotatedMean, KERSIZE * KERSIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);
                        alreadyBcasted[4] = 1;
                    }
                    break; 
            case 5: if (alreadyBcasted[5] == 0) {
                        MPI_Bcast(rotatedEmboss, KERSIZE * KERSIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);
                        alreadyBcasted[5] = 1;
                    }
                    break; 
        }
        i++;
    }

    if (rank == MASTER) { 
        // Send no of rows to be processed and offsets to each other proc
        averageRows = height / nProcesses;
        restRows = height % nProcesses;
        offset[1] = averageRows + 1;

        for (int i = 1; i < nProcesses; i++) {
            rowsPerProc[i] = (restRows >= i) ? averageRows + 1 : averageRows;
            MPI_Send(&offset[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&rowsPerProc[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            
            offset[i + 1] = offset[i] + rowsPerProc[i];
        }

        // Save no of rows to be processed in MASTER proc
        rowsPerProcLoc = averageRows;

        // For each filter application, send input image matrix
        for (int k = 0; k < filtSize; k++) {
            for (int i = 1; i < nProcesses; i++) {
                // Send matrix according to image type
                if (strcmp(type, BWIMG) == 0) {
                    if (isT2forResult) {
                        // Send 2 more rows each time, for above-below neighbours
                        MPI_Send(&bwMatrix1[offset[i] - 1][0], (rowsPerProc[i] + 2) * (width + 2), 
                            MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
                    } else {
                        MPI_Send(&bwMatrix2[offset[i] - 1][0], (rowsPerProc[i] + 2) * (width + 2), 
                            MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
                    }
                } else if (strcmp(type, COLORIMG) == 0) {
                    if (isT2forResult) {
                        MPI_Send(&rgbMatrix1[offset[i] - 1][0], (rowsPerProc[i] + 2) * (width + 2), 
                            mpi_rgbmatrix, i, 0, MPI_COMM_WORLD);
                    } else {
                        MPI_Send(&rgbMatrix2[offset[i] - 1][0], (rowsPerProc[i] + 2) * (width + 2), 
                            mpi_rgbmatrix, i, 0, MPI_COMM_WORLD);
                    }
                }
            }
            
            // Apply filter for rows processed by MASTER proc
            if (strcmp(type, BWIMG) == 0) {
                if (isT2forResult) {
                    // Apply according to image type and I/O matrices
                    switch (filtToApply[k]) { 
                        case 1: applyFilter(1, width, rotatedSmooth, bwMatrix1, bwMatrix2, 
                                    NULL, NULL, rowsPerProcLoc);
                                break; 
                        case 2: applyFilter(1, width, rotatedBlur, bwMatrix1, bwMatrix2, 
                                    NULL, NULL, rowsPerProcLoc);
                                break; 
                        case 3: applyFilter(1, width, rotatedSharpen, bwMatrix1, bwMatrix2, 
                                    NULL, NULL, rowsPerProcLoc);
                                break;
                        case 4: applyFilter(1, width, rotatedMean, bwMatrix1, bwMatrix2, 
                                    NULL, NULL, rowsPerProcLoc);
                                break;
                        case 5: applyFilter(1, width, rotatedEmboss, bwMatrix1, bwMatrix2, 
                                    NULL, NULL, rowsPerProcLoc);
                                break;
                    }
                } else {
                    switch (filtToApply[k]) { 
                        case 1: applyFilter(1, width, rotatedSmooth, bwMatrix2, bwMatrix1, 
                                    NULL, NULL, rowsPerProcLoc);
                                break; 
                        case 2: applyFilter(1, width, rotatedBlur, bwMatrix2, bwMatrix1, 
                                    NULL, NULL, rowsPerProcLoc);
                                break; 
                        case 3: applyFilter(1, width, rotatedSharpen, bwMatrix2, bwMatrix1, 
                                    NULL, NULL, rowsPerProcLoc);
                                break;
                        case 4: applyFilter(1, width, rotatedMean, bwMatrix2, bwMatrix1, 
                                    NULL, NULL, rowsPerProcLoc);
                                break;
                        case 5: applyFilter(1, width, rotatedEmboss, bwMatrix2, bwMatrix1, 
                                    NULL, NULL, rowsPerProcLoc);
                                break;
                    }
                }
            } else if (strcmp(type, COLORIMG) == 0) {
                if (isT2forResult) {
                    switch (filtToApply[k]) { 
                        case 1: applyFilter(1, width, rotatedSmooth, NULL, NULL, rgbMatrix1, 
                                    rgbMatrix2, rowsPerProcLoc);
                                break; 
                        case 2: applyFilter(1, width, rotatedBlur, NULL, NULL, rgbMatrix1, 
                                    rgbMatrix2, rowsPerProcLoc);
                                break; 
                        case 3: applyFilter(1, width, rotatedSharpen, NULL, NULL, rgbMatrix1, 
                                    rgbMatrix2, rowsPerProcLoc);
                                break;
                        case 4: applyFilter(1, width, rotatedMean, NULL, NULL, rgbMatrix1, 
                                    rgbMatrix2, rowsPerProcLoc);
                                break;
                        case 5: applyFilter(1, width, rotatedEmboss, NULL, NULL, rgbMatrix1, 
                                    rgbMatrix2, rowsPerProcLoc);
                                break;
                    }
                } else {
                    switch (filtToApply[k]) { 
                        case 1: applyFilter(1, width, rotatedSmooth, NULL, NULL, rgbMatrix2, 
                                    rgbMatrix1, rowsPerProcLoc);
                                break; 
                        case 2: applyFilter(1, width, rotatedBlur, NULL, NULL, rgbMatrix2, 
                                    rgbMatrix1, rowsPerProcLoc);
                                break; 
                        case 3: applyFilter(1, width, rotatedSharpen, NULL, NULL, rgbMatrix2, 
                                    rgbMatrix1, rowsPerProcLoc);
                                break;
                        case 4: applyFilter(1, width, rotatedMean, NULL, NULL, rgbMatrix2, 
                                    rgbMatrix1, rowsPerProcLoc);
                                break;
                        case 5: applyFilter(1, width, rotatedEmboss, NULL, NULL, rgbMatrix2, 
                                    rgbMatrix1, rowsPerProcLoc);
                                break;
                    }
                }
            }

            // Receive processed rows from other procs
            for (int i = 1; i < nProcesses; i++) {
                if (strcmp(type, BWIMG) == 0) {
                    if (isT2forResult) {
                        // Store according to output matrix for this step
                        MPI_Recv(&bwMatrix2[offset[i]][0], rowsPerProc[i] * (width + 2), MPI_UNSIGNED_CHAR, 
                            i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    } else {
                        MPI_Recv(&bwMatrix1[offset[i]][0], rowsPerProc[i] * (width + 2), MPI_UNSIGNED_CHAR, 
                            i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                } else if (strcmp(type, COLORIMG) == 0) {
                    if (isT2forResult) {
                        MPI_Recv(&rgbMatrix2[offset[i]][0], rowsPerProc[i] * (width + 2), mpi_rgbmatrix, 
                            i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    } else {
                        MPI_Recv(&rgbMatrix1[offset[i]][0], rowsPerProc[i] * (width + 2), mpi_rgbmatrix, 
                            i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }

            // Change output to input matrix for next stept
            if (isT2forResult) {
                isT2forResult = 0;
            } else {
                isT2forResult = 1;
            }
            
        }

        // Save final image to output file
        if (isT2forResult) {
            saveFile(file_out, width, height, rgbMatrix1, bwMatrix1, vmax, type); 
        } else {
            saveFile(file_out, width, height, rgbMatrix2, bwMatrix2, vmax, type); 
        }

    } else {

        // Receive number of rows to be processed and offset
        MPI_Recv(&offsetLoc, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rowsPerProcLoc, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Receive input rows from MASTER
        for (int k = 0; k < filtSize; k++) {
            if (strcmp(type, BWIMG) == 0) {
                MPI_Recv(&bwMatrix1[offsetLoc - 1][0], (rowsPerProcLoc + 2) * (width + 2), MPI_UNSIGNED_CHAR, 
                    MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else if (strcmp(type, COLORIMG) == 0) {
                MPI_Recv(&rgbMatrix1[offsetLoc - 1][0], (rowsPerProcLoc + 2) * (width + 2), mpi_rgbmatrix, 
                    MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // Apply filter according to image type
            if (strcmp(type, BWIMG) == 0) {
                switch (filtToApply[k]) { 
                    case 1: applyFilter(offsetLoc, width, rotatedSmooth, bwMatrix1, bwMatrix2, 
                                NULL, NULL, rowsPerProcLoc);
                            break; 
                    case 2: applyFilter(offsetLoc, width, rotatedBlur, bwMatrix1, bwMatrix2, 
                                NULL, NULL, rowsPerProcLoc);
                            break; 
                    case 3: applyFilter(offsetLoc, width, rotatedSharpen, bwMatrix1, bwMatrix2, 
                                NULL, NULL, rowsPerProcLoc);
                            break;
                    case 4: applyFilter(offsetLoc, width, rotatedMean, bwMatrix1, bwMatrix2, 
                                NULL, NULL, rowsPerProcLoc);
                            break;
                    case 5: applyFilter(offsetLoc, width, rotatedEmboss, bwMatrix1, bwMatrix2, 
                                NULL, NULL, rowsPerProcLoc);
                            break;
                }
            } else if (strcmp(type, COLORIMG) == 0) {
                switch (filtToApply[k]) { 
                    case 1: applyFilter(offsetLoc, width, rotatedSmooth, NULL, NULL, rgbMatrix1, 
                                rgbMatrix2, rowsPerProcLoc);
                            break; 
                    case 2: applyFilter(offsetLoc, width, rotatedBlur, NULL, NULL, rgbMatrix1, 
                                rgbMatrix2, rowsPerProcLoc);
                            break; 
                    case 3: applyFilter(offsetLoc, width, rotatedSharpen, NULL, NULL, rgbMatrix1, 
                                rgbMatrix2, rowsPerProcLoc);
                            break;
                    case 4: applyFilter(offsetLoc, width, rotatedMean, NULL, NULL, rgbMatrix1, 
                                rgbMatrix2, rowsPerProcLoc);
                            break;
                    case 5: applyFilter(offsetLoc, width, rotatedEmboss, NULL, NULL, rgbMatrix1, 
                                rgbMatrix2, rowsPerProcLoc);
                            break;
                }
            }

            // Send processed rows back to MASTER proc
            if (strcmp(type, BWIMG) == 0) {
                MPI_Send(&bwMatrix2[offsetLoc][0], rowsPerProcLoc * (width + 2), MPI_UNSIGNED_CHAR, MASTER, 
                    0, MPI_COMM_WORLD);
            } else if (strcmp(type, COLORIMG) == 0) {
                MPI_Send(&rgbMatrix2[offsetLoc][0], rowsPerProcLoc * (width + 2), mpi_rgbmatrix, MASTER, 
                    0, MPI_COMM_WORLD);
            }
        }
    }

    // Free dynamically allocated memory 
    if (strcmp(type, BWIMG) == 0) {
        free(bwMatrix1[0]);
        free(bwMatrix1);

        free(bwMatrix2[0]);
        free(bwMatrix2);
    } else if (strcmp(type, COLORIMG) == 0) {
        free(rgbMatrix1[0]);
        free(rgbMatrix1);

        free(rgbMatrix2[0]);
        free(rgbMatrix2);
    }
    MPI_Type_free(&mpi_rgbmatrix);

    MPI_Finalize();
    return 0;
}