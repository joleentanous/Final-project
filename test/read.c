#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INITIAL_BUFFER_SIZE 128

char* read_line(FILE* file) {
    size_t buffer_size = INITIAL_BUFFER_SIZE;
    size_t length = 0;
    char* buffer = (char*)malloc(buffer_size * sizeof(char));
    int ch;

    if (buffer == NULL) {
        perror("Unable to allocate buffer");
        exit(EXIT_FAILURE);
    }

    while ((ch = fgetc(file)) != '\n' && ch != EOF) {
        buffer[length++] = ch;
        if (length == buffer_size) {
            buffer_size *= 2;
            buffer = (char*)realloc(buffer, buffer_size * sizeof(char));
            if (buffer == NULL) {
                perror("Unable to reallocate buffer");
                exit(EXIT_FAILURE);
            }
        }
    }

    if (length == 0 && ch == EOF) {
        free(buffer);
        return NULL;
    }

    buffer[length] = '\0';
    return buffer;
}

void count_rows_and_columns(const char* file_name, int* rows, int* cols) {
    FILE *file;
    char *line;
    int temp_cols;
    char *token;

    file = fopen(file_name, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    *rows = 0;
    *cols = 0;

    while ((line = read_line(file)) != NULL) {
        (*rows)++;
        token = strtok(line, ",");
        temp_cols = 0;
        while (token != NULL) {
            temp_cols++;
            token = strtok(NULL, ",");
        }
        if (*cols == 0) {
            *cols = temp_cols;
        } else if (temp_cols != *cols) {
            fprintf(stderr, "Inconsistent number of columns\n");
            free(line);
            exit(EXIT_FAILURE);
        }
        free(line);
    }

    fclose(file);
}

void read_data(const char* file_name, double* data, int rows, int cols) {
    FILE *file;
    char *line;
    int i = 0;
    char *token;
    rows = rows;
    cols = cols;
    file = fopen(file_name, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    while ((line = read_line(file)) != NULL) {
        token = strtok(line, ",");
        while (token != NULL) {
            data[i++] = atof(token);
            token = strtok(NULL, ",");
        }
        free(line);
    }

    fclose(file);
}

int main(int argc, char *argv[]) {
    const char* operation;
    const char* file_name;
    int N, d;
    double* data;
    int i;

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <operation> <file_name>\n", argv[0]);
        return EXIT_FAILURE;
    }

    operation = argv[1];
    file_name = argv[2];

    count_rows_and_columns(file_name, &N, &d);

    data = (double*)malloc((size_t)(N * d) * sizeof(double));
    if (data == NULL) {
        perror("Memory allocation failed");
        return EXIT_FAILURE;
    }

    read_data(file_name, data, N, d);

    printf("N = %d, d = %d\n", N, d);
    printf("Data:\n");
    for (i = 0; i < N * d; i++) {
        printf("%f ", data[i]);
        if ((i + 1) % d == 0) {
            printf("\n");
        }
    }
    printf("%s\n", operation);
    free(data);
    return EXIT_SUCCESS;
}
