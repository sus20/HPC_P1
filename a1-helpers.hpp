/**
 *   This file contains a set of helper function. 
 *   Comment: You probably do need to change any of these.
*/
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

double** allocate(int height, int width, const double& val = 0) {    
    double** ptr = new double*[height]; 
    double* mem = new double[height*width]{ val }; 

    for (unsigned i = 0; i < height; ++i, mem += width)
        ptr[i] = mem;

    return ptr;
}

void deallocate(double** data) {
   delete [] data[0];  
   delete [] data;     
}

/**
 * Matrix data structure
 * 2D => (y-dim, x-dim)
 * Init:
 *   Mat mat(800, 600);

 * Access element at position (34, 53):
 *  mat(34,53)
 * 
 * Access height and width with: 
 * mat.height and mat.width
 * 
*/

struct Mat {
    double** data = nullptr;
    
    int height;
    int width;
    
    
    Mat (int height, int width, const double& val = 0) 
        : height(height), width(width), data(nullptr)
    {   
        if ( height > 0 && width >  0)
            data = allocate(height, width, 0);
    }

    ~Mat() {
        if (data) {
            delete [] data[0];  
            delete [] data; 
        }
    }
    // Mat() : height(0), width(0), data(0) {}

    // Mat(const Mat& rhs) : height(rhs.height), width(rhs.width), data(height*width, val) 
    // {
    //     std::copy(rhs.data().data[0][0], &rhs.data[height-1][width], &data[0][0]);
    // }

    double& operator()(int h, int w) {
        return data[h][w];
    }

    void swap(Mat& right)
    {
        std::swap(this->data, right.data);
        std::swap(this->width, right.width);
        std::swap(this->height, right.height);
    }

    double* operator[](unsigned row)
    {
        return data[row];
    }

    const double* operator[](unsigned row) const
    {
        return data[row];
    }
    
    std::string to_string() {
        std::stringstream ss;
        for (int j = 0; j < height; ++j) {
            for (int i = 0; i < width; ++i) {
                ss << std::fixed << std::setfill(' ') << std::right << std::setw(11) << std::setprecision(6)  << this->operator()(j,i) << " ";
                
            }
            ss << std::endl;
        }
        
        return ss.str();
    }
};

/**
 * Process command line arguments
 * Note: N and M are required arguments
 * 
*/
void process_input(int argc, char **argv, int& N, int& M, int& nx, int& ny, int& max_iterations, double& epsilon) {
    int input_check = 0;

    // process input arguments
    for (int i = 0; i < argc; ++i) {
        if ( std::string(argv[i]).compare("--n") == 0 ) {
            N=atoi(argv[++i]);
            input_check++;
        }
        if ( std::string(argv[i]).compare("--m") == 0 ) {
            M=atoi(argv[++i]);
            input_check++;
        }
        if ( std::string(argv[i]).compare("--nx") == 0 ) {
            nx=atoi(argv[++i]);
            input_check++;
        }
        if ( std::string(argv[i]).compare("--ny") == 0 ) {
            ny=atoi(argv[++i]);
            input_check++;
        }
        if ( std::string(argv[i]).compare("--max-iterations") == 0 ) {
            max_iterations=atoi(argv[++i]);
        }
        if ( std::string(argv[i]).compare("--epsilon") == 0 ) {
            epsilon=atof(argv[++i]);
        }
    }
    if ( input_check < 2 ) {
        std::cout << "Usage: --n <columns> --m <rows> [--max-iterations <int>]" << std::endl;
        exit(-1);
    }
}

/**
 * Print given Matrix of type Mat to std out 
 * 
*/
void print_mat(Mat& m) {
    std::stringstream ss;
    for (int j = 0; j < m.height; ++j) {
        for (int i = 0; i < m.width; ++i) {
            ss << m(j,i) << " ";
            
        }
        ss << std::endl;
    }
    std::cout << ss.str();
}

/**
 * Checks if all values of two input matrices are equal.
 * 
 * @param[in] a
 * @param[in] b
 * 
 * Return true if two matrices of type Mat have the same values.
*/
bool verify(Mat& a, Mat& b) {
    if ( a.height != b.height || a.width != b.width ) 
        return false;
    
    for (int j=0;j<a.height;++j) {
        for (int i=0;i<a.width;++i) {
            if ( std::abs( a(j,i) -  b(j,i)) > std::pow(10, -8) ) {
                return false; 
            }
        }
    }

    return true;
}

/**
 * Jacobi iterative solver for Heat Equation - sequential version
 * @param[in] max_iterations
 * @param[in] epsilon
 * @param[inout] U
 * @param[inout] iteration_count
*/
void heat2d_sequential(int max_iterations, double epsilon, Mat& U, int& iteration_count) {
    int i, j;
    double diffnorm;

    int M = U.height, N = U.width;

    // allocate another 2D array
    Mat W(M,N); 

    // Init & Boundary
    for (j = 0; j < M; ++j) {
        for (i = 0; i < N; ++i) {
            W[j][i] = U[j][i] = 0.0;
        }

        W[j][0] = U[j][0] = 20.0;
        W[j][N-1] = U[j][N-1] = 40.0;
    }

    for (i = 0; i < N; ++i) {
        W[M - 1][i] = U[M - 1][i] = 100.0;
    }
    

    // End init

    iteration_count = 0;
    do
    {
        iteration_count++;
        diffnorm = 0.0;
        
        /* Compute new values (but not on boundary) */
        for (i = 1; i < M - 1; ++i){
            for (j = 1; j < N - 1; ++j)
            {
                W(i,j) = (U(i,j + 1) + U(i,j - 1) + U(i + 1,j) + U(i - 1,j)) * 0.25;
                diffnorm += (W(i,j) - U(i,j))*(W(i,j) - U(i,j));
            }
        }

        // Only transfer the interior points
        for (i = 1; i < M - 1; ++i)
            for (j = 1; j < N - 1; ++j)
                U(i,j) = W(i,j);

        diffnorm = sqrt(diffnorm);
    } while (epsilon <= diffnorm  && iteration_count < max_iterations);
}

/**
 * Save given Matrix of type Mat to disk in PPM image format
 * Note: N and M are required arguments
 * 
*/
void save_to_disk(Mat& U, std::string filename) {
    std::ofstream ofs(filename, std::ofstream::out);
    ofs << U.width << "\n"
        << U.height << "\n";
    
    for (int i = 0; i < U.height; ++i)
    {
        for (int j = 0; j < U.width; ++j)
        {
            ofs << std::fixed << std::setfill(' ') << std::right << std::setw(11) << std::setprecision(6) <<  U(i,j) << " ";
        }
        // ofs << "\n";
    }
    ofs.close();
}
