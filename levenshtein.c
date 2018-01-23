//
// Created by nicolet on 19/01/18.
//


int levenshtein(char* X, char *Y, m){
    int **D;
    int i,j;

    D = malloc(sizeof(D)*(n + 1));
    for(i = 0; i< n + 1; i++){
        D[i] = malloc(sizeof(D[i])*(n + 1));
    }

    for(i = 0; i < n+1; i++){
        D[0][i] = i;
        D[i][0] = i;
    }
#pragma scop
    for(i = 1; i < n + 1; i++){
        for(j = 1; j < n + 1; j++){
            offs = X[i-1] == Y[j-1] ? 2 : 0;
            D[i][j] = D[i-1][j] + D[i][j-1] + D[i-1][j-1] + offs;
        }
    }
#pragma endscop
    int res = D[n][n];

    for(i = 0; i< n; i++){
        free(D[i]);
    }
    free(D);

    return res;
}

