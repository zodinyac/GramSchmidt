#include <iostream>
#include <math.h>
#include <mpi.h>
#include <map>
#include <queue>
#include <fstream>
#include <iomanip>

using namespace std;

const double eps = 1e-7;
enum {
    NO_MORE_TASKS = -1,
    GET_TASK_NUMBER,
    GET_VECTOR_STATUS,
    SEND_VECTOR_STATUS,
    IS_ZERO_VECTOR,
    IS_NON_ZERO_VECTOR
};

enum {
    NULL_VECTOR = -2,
    ZERO_VECTOR = -1
};

void assignment(double *destination, double *source, int length) {
    while (length > 0) {
        *destination = *source;
        ++source;
        ++destination;
        --length;
    }
}

void proj_and_minus(double *destination, double *a, double *b, int length) {
    double top(0), bottom(0);
    double *tmp_b = b;
    for (int i = 0;  i < length; ++i, ++a, ++tmp_b) {
        top += *a * *tmp_b;
        bottom += *tmp_b * *tmp_b;
    }
    double s = top / bottom;
    for (int i = 0;  i < length; ++i, ++destination, ++b)
        *destination -= s * *b;
}

bool is_non_zero(double *a, int length) {
    for (int i = 0; i < length; ++i)
        if (fabs(a[i]) > eps)
            return true;
    return false;
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        cout << "Usage: gramschmidt <input.dat> <output.dat>." << endl;
        return 1;
    }

    MPI_Status status;
    int ret, rank, size;
    double algo_time, write_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_File input, output;
    ret = MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &input);
    if (ret != MPI_SUCCESS) {
        if (rank == 0) {
            printf("File '%s' can't be opened.\n", argv[1]);
        }
        MPI_Finalize();
        return 2;
    }
    ret = MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output);
    if (ret != MPI_SUCCESS) {
        if (rank == 0) {
            printf("File '%s' can't be opened.\n", argv[1]);
        }
        MPI_Finalize();
        return 2;
    }

    int m, n;
    MPI_File_read(input, &m, 1, MPI_INT, &status);
    MPI_File_read(input, &n, 1, MPI_INT, &status);

    if (!rank) {
        cout << "Proc count: " << size << endl;
        cout << "Matrix size: " << m << " x " << n << endl;
    }
    algo_time = MPI_Wtime();

    if (size == 1) {
        // serial version
        double **a, **b;
        bool *non_zero_b;
        a = new double*[m];
        b = new double*[m];
        non_zero_b = new bool[m];
        for (int i = 0; i < m; ++i) {
            a[i] = new double[n];
            b[i] = new double[n];
        }

        for (int i = 0; i < m; ++i)
            MPI_File_read(input, a[i], n, MPI_DOUBLE, &status);

        assignment(b[0], a[0], n);
        non_zero_b[0] = is_non_zero(b[0], n);
        for (int i = 1; i < m; ++i) {
            assignment(b[i], a[i], n);
            for (int j = 0; j < i; ++j)
                if (non_zero_b[j])
                    //proj_and_minus(b[i], a[i], b[j], n); // default scheme
                    proj_and_minus(b[i], b[i], b[j], n); // pro scheme
            non_zero_b[i] = is_non_zero(b[i], n);
        }

        algo_time = MPI_Wtime() - algo_time;
        write_time = MPI_Wtime();

        int b_m = 0;
        for (int i = 0; i < m; ++i)
            b_m += non_zero_b[i];
        MPI_File_write(output, &b_m, 1, MPI_INT, &status);
        MPI_File_write(output, &n, 1, MPI_INT, &status);
        for (int i = 0; i < m; ++i)
            if (non_zero_b[i])
                MPI_File_write(output, b[i], n, MPI_DOUBLE, &status);

        write_time = MPI_Wtime() - write_time;
        cout << fixed << setprecision((int)round(log(1 / MPI_Wtick()) / log(10)));
        cout << "Algo time: " << algo_time << endl;
        cout << "Write time: " << write_time << endl;

        for (int i = 0; i < m; ++i) {
            delete[] a[i];
            delete[] b[i];
        }
        delete[] a;
        delete[] b;
        delete[] non_zero_b;
        MPI_File_close(&input);
        MPI_File_close(&output);
        MPI_Finalize();
        return 0;
    } else {
        // parallel version
        double **b = new double*[m];
        int b_m = 0;
        int *b_ind2ind = new int[m];
        for (int i = 0; i < m; ++i)
            b_ind2ind[i] = NULL_VECTOR;
        if (!rank) {
            //master
            int i = 0, command, ended_procs = 0, ind;
            map<int, queue<int> > waiting_procs;
            while (true) {
                MPI_Recv(&command, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                switch (command) {
                case GET_TASK_NUMBER:
                    if (i < n) {
                        MPI_Send(&command, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                        MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                        ++i;
                    } else {
                        command = NO_MORE_TASKS;
                        MPI_Send(&command, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                        ++ended_procs;
                        if (ended_procs == size - 1) {
                            // exit
                            algo_time = MPI_Wtime() - algo_time;
                            write_time = MPI_Wtime();

                            MPI_File_write(output, &b_m, 1, MPI_INT, &status);
                            MPI_File_write(output, &n, 1, MPI_INT, &status);
                            for (int i = 0; i < b_m; ++i) {
                                MPI_File_write(output, b[i], n, MPI_DOUBLE, &status);
                                delete[] b[i];
                            }

                            write_time = MPI_Wtime() - write_time;
                            cout << fixed << setprecision((int)round(log(1 / MPI_Wtick()) / log(10)));
                            cout << "Algo time: " << algo_time << endl;
                            cout << "Write time: " << write_time << endl;

                            delete[] b;
                            delete[] b_ind2ind;
                            MPI_File_close(&input);
                            MPI_File_close(&output);
                            MPI_Finalize();
                            return 0;
                        }
                    }
                    break;
                case GET_VECTOR_STATUS:
                    MPI_Recv(&ind, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    switch (b_ind2ind[ind]) {
                    case NULL_VECTOR:
                        waiting_procs[ind].push(status.MPI_SOURCE);
                        break;
                    case ZERO_VECTOR:
                        command = IS_ZERO_VECTOR;
                        MPI_Send(&command, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                        break;
                    default:
                        command = IS_NON_ZERO_VECTOR;
                        MPI_Send(&command, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                        MPI_Send(b[b_ind2ind[ind]], n, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                    }
                    break;
                case SEND_VECTOR_STATUS:
                    MPI_Recv(&ind, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    MPI_Recv(&command, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    if (command == IS_NON_ZERO_VECTOR) {
                        b[b_m] = new double[n];
                        MPI_Recv(b[b_m], n, MPI_DOUBLE, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        b_ind2ind[ind] = b_m;
                        while (!(waiting_procs[ind].empty())) {
                            MPI_Send(&command, 1, MPI_INT, waiting_procs[ind].front(), 0, MPI_COMM_WORLD);
                            MPI_Send(b[b_m], n, MPI_DOUBLE, waiting_procs[ind].front(), 0, MPI_COMM_WORLD);
                            waiting_procs[ind].pop();
                        }
                        ++b_m;
                    } else {
                        b_ind2ind[ind] = ZERO_VECTOR;
                        while (!(waiting_procs[ind].empty())) {
                            MPI_Send(&command, 1, MPI_INT, waiting_procs[ind].front(), 0, MPI_COMM_WORLD);
                            waiting_procs[ind].pop();
                        }
                    }
                }
            }
        } else {
            //slaves
            double *b_tmp;
            int command, task;
            while (true) {
                command = GET_TASK_NUMBER;
                MPI_Sendrecv_replace(&command, 1, MPI_INT, 0, 0, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                switch (command) {
                case NO_MORE_TASKS:
                    // exit
                    for (int i = 0; i < b_m; ++i)
                        delete[] b[i];
                    delete[] b;
                    delete[] b_ind2ind;
                    MPI_File_close(&input);
                    MPI_File_close(&output);
                    MPI_Finalize();
                    return 0;
                case GET_TASK_NUMBER:
                    MPI_Recv(&task, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    MPI_File_seek(input, 2 * sizeof(int) + n * task * sizeof(double), MPI_SEEK_SET);
                    b_tmp = new double[n];
                    MPI_File_read(input, b_tmp, n, MPI_DOUBLE, &status); // bi = ai
                    for (int j = 0; j < task; ++j) {
                        if (b_ind2ind[j] == NULL_VECTOR) {
                            command = GET_VECTOR_STATUS;
                            MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                            command = j;
                            MPI_Sendrecv_replace(&command, 1, MPI_INT, 0, 0, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                            if (command == IS_NON_ZERO_VECTOR) {
                                b[b_m] = new double[n];
                                MPI_Recv(b[b_m], n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // tmp = bj
                                b_ind2ind[j] = b_m;
                                ++b_m;
                                proj_and_minus(b_tmp, b_tmp, b[b_ind2ind[j]], n); // pro scheme
                            }
                        } else if (b_ind2ind[j] != ZERO_VECTOR)
                            proj_and_minus(b_tmp, b_tmp, b[b_ind2ind[j]], n); // pro scheme
                    }
                    command = SEND_VECTOR_STATUS;
                    MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    command = task;
                    MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    if (is_non_zero(b_tmp, n)) {
                        command = IS_NON_ZERO_VECTOR;
                        MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                        MPI_Send(b_tmp, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                        b[b_m] = b_tmp;
                        b_ind2ind[task] = b_m;
                        ++b_m;
                    } else {
                        command = IS_ZERO_VECTOR;
                        MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                        b_ind2ind[task] = ZERO_VECTOR;
                        delete[] b_tmp;
                    }
                    break;
                }
            }
        }
    }
}
