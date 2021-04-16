// Templated functions to perform convolution operations on type T vectors
// Either std::vector<std::vector<gray>> or std::vector<std::vector<ppm_pixel>>
// Depends on what std::holds_alternative evaluates to.
#ifndef UTILS_H
#define UTILS_H

#include <sstream>
#include <string>


template<typename T>
double ninety_degree_flip(T& image, const unsigned short n_threads) {
  

  // Much like matrix transpose, right?
  T convoluted_vec(image[0].size());
  // TODO: May need to change the dimensions of this vector
  for (auto& internal_vec : convoluted_vec) {
    internal_vec.resize(image.size());
  }

  double start_parallel(omp_get_wtime());
  // TODO: Define parallel region
  #pragma omp parallel for num_threads(n_threads)
  for (int i = 0; i < image.size(); ++i) {
    for (int j = 0; j < image[0].size(); ++j) {
      convoluted_vec[j][(image.size()-1)-i] = image[i][j]; 
    }
  }
 

  // end of parllel region
  double end_parallel(omp_get_wtime());
  double total_parallel(end_parallel - start_parallel);
  // moving convoluted_vec into image
  image = std::move(convoluted_vec);
  return  total_parallel;
}

template<typename T>
double vertical_flip(T& image, const unsigned short n_threads) {

  double start_parallel(omp_get_wtime());
  // TODO: Define parallel region
  #pragma omp parallel for num_threads(n_threads)
  for (int j = 0; j < image[0].size(); ++j) {
    for (int i = 0; i < (image.size() / 2); ++i) {
      auto temp = image[i][j];
      image[i][j] = image[(image.size()-1)-i][j];
      image[(image.size()-1)-i][j] = temp;
    }
  } 
  // end of parllel region
  double end_parallel(omp_get_wtime());
  double total_parallel(end_parallel - start_parallel);
  return total_parallel;
}


template<typename T>
double horizontal_flip(T& image, const unsigned short n_threads) {

  double start_parallel(omp_get_wtime());
  #pragma omp parallel for num_threads(n_threads) 
  for (int i = 0; i < image.size(); ++i) {
    for (int j = 0; j < (image[0].size() / 2); ++j) {
      auto temp = image[i][j];
      image[i][j] = image[i][(image[i].size()-1)-j];
      image[i][(image[i].size()-1)-j] = temp;
    }
  }
  // end of parllel region
  double end_parallel(omp_get_wtime());
  double total_parallel(end_parallel - start_parallel);
  return total_parallel;
}

template<typename T>
void print_info(char** argv, T& image, const unsigned short n_threads, const double parallel_time, const double total_time) {
  std::stringstream ss;
  std::string program_name(argv[0]);
  program_name = program_name.substr(program_name.find_last_of('/')+1);
  ss << program_name << ',' << image[0].size() << ',' << image.size() << ',' <<
    n_threads << ',' << parallel_time << ',' << total_time;
  std::cout << ss.str() << std::endl;
}


#endif 
