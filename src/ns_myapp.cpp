// myapp.cpp
#include <stdexcept>
#include <sstream>
#include <thread>
#include <random>
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include <system_error>
#include <random>
#include <omp.h>
#include "ns_seed.hpp"   // seed::get_global_seed()

//#####################################################################################
//#####################################################################################
// @namespace myapp
// @brief namespace for miscellaneous functions
//#####################################################################################
//#####################################################################################


// open file and get fout with safety
/// @hidecallgraph
/// @hidecallergraph
FILE* myapp::get_fout( const std::filesystem::path& pathout )
{
  if(pathout.empty()) THROW_ERROR("fname is not assigned");
  FILE *fout = fopen(pathout.c_str(),"wt");
  if(fout == NULL) THROW_ERROR("Cannot open file : " + pathout.string());
  return fout;
}

FILE* myapp::get_fout_binary( const std::filesystem::path &filepath )
{
  if(filepath.empty()) THROW_ERROR("fname is not assigned");
  FILE *fout = fopen(filepath.c_str(), "wb");
  if(fout == NULL) THROW_ERROR("Cannot open file : " + filepath.string());
  return fout;
}

// @brief close file with safety
void myapp::close( FILE* fout, const std::filesystem::path& pathout )
{
  if(fclose(fout) == EOF){
    THROW_ERROR("myapp::close: fclose failed. fname={}", pathout.string());
  }
}

// @brief check file exists or not
void myapp::filecheck(const std::filesystem::path& pathout )
{
  if( std::filesystem::exists(pathout) == false ){
    THROW_ERROR("myapp::filecheck: File does not exist. path={}", pathout.string());
  }else{
    LOG_INFO("File exists. path={}", pathout.string());
  }
}


// 2022-10-31 17:08:41
// https://qiita.com/iseki-masaya/items/70b4ee6e0877d12dafa8
std::vector< std::string > myapp::split_naive(const std::string &s, const char char_delim)
{
  std::vector< std::string > elems;
  std::string item;
  for (char ch: s) {
    if (ch == char_delim) {
      if ( !item.empty() ) elems.push_back(item);
      item.clear();
    }
    else {
      item += ch;
    }
  }
  if ( !item.empty() ) elems.push_back(item);
  return elems;
}

// Avoid returning empty strings even with multiple spaces
std::vector<std::string> myapp::split(
  const std::string &s, const char char_delim, const char char_commentout)
{
  std::vector<std::string> elems;
  std::string item;
  bool isCommentOut = false;

  for (char ch : s) {
    if (ch == char_commentout) {
      isCommentOut = true;
    }

    if (!isCommentOut) {
      if (ch == char_delim) {
        if (!item.empty()) {
          elems.push_back(item);
          item.clear();
        }
      }
      else {
        item += ch;
      }
    }
  }
  if (!isCommentOut && !item.empty()) elems.push_back(item);
  return elems;
}

// 2023-04-10 15:21:54
// Round to n decimal places
float myapp::round_n( float number, const double n)
{
  if( n < 1 )
    THROW_ERROR("float myapp::round_n ,  n < 1 ");
  number = number * pow(10,n-1); // Multiply value to be rounded by 10^(n-1).
  number = round(number); // Round to nearest integer.
  number /= pow(10, n-1); // Divide by 10^(n-1).
  return number;
}

double myapp::round_n( double number, const double n)
{
  if( n < 1 )
    THROW_ERROR("double myapp::round_n ,  n < 1 ");
  number = number * pow(10,n-1); // Multiply value to be rounded by 10^(n-1).
  number = round(number); // Round to nearest integer.
  number /= pow(10, n-1); // Divide by 10^(n-1).
  return number;
}

namespace {
// Reproducible per-thread RNG (global_seed + thread_id)
std::mt19937 make_poisson_gen()
{
  std::seed_seq seq{seed::get_global_seed(), static_cast<unsigned>(omp_get_thread_num())};
  return std::mt19937(seq);
}
} // anonymous

// Get random value from Poisson distribution
int myapp::poisson(const double mean)
{
  if (mean < 0)
    THROW_ERROR("myapp::poisson_distribution, mean<0");
  thread_local std::mt19937 gen = make_poisson_gen();
  std::poisson_distribution<int> dist(mean);
  return dist(gen);
}

// Count digits
// https://www.delftstack.com/ja/howto/cpp/number-of-digits-in-a-number-cpp/
// size_t myapp::count_digits(const double number)
// {
//   std::string tmp;
//   tmp = std::to_string(number);
//   tmp.erase(std::remove_if(tmp.begin(), tmp.end(), ::ispunct), tmp.end());
//   return tmp.size();
// }

// by  chatGPT4
// 2023-04-10 16:43:06
// arranged by me
// Count decimal places
// The smaller seido_keta is, the more accurate the result.
int myapp::countDecimalPlaces(double number, const int seido_keta)
{
  const double factor = pow(10,seido_keta);
  const double small = factor*fabs(number)*std::numeric_limits<double>::epsilon(); // Set appropriate precision
  int decimalPlaces = 0;
  while (std::abs(number - static_cast<long long>(number)) > small) {
      number *= 10;
      decimalPlaces++;
  }
  return decimalPlaces;
}

double myapp::log10(const double x)
{
  if(x<0)
    THROW_ERROR("myapp::log10 , x< 0");
  return ::log10(x); // To call log10 not in myapp:: namespace, use ::log10().
}

/// @brief Convert number to string with specified digit format
std::string myapp::formatNumber(int digits, int number) {
  char buffer[100]; // Buffer of sufficient size
  std::string format = "%0" + std::to_string(digits) + "d";
  sprintf(buffer, format.c_str(), number); // Format number and write to string
  return std::string(buffer); // Convert to std::string and return
}

// 2022-08-22 11:38:48
// bool_zero or not for double
// factor is for std::numeric_limits<double>::epsilon()
bool myapp::is_zero( const double value, const double factor )
{
  const double small = std::numeric_limits<double>::epsilon()*fabs(factor);
  if( fabs(value) < small ) return true;
  return false;
}

// 2024-09-20 16:16:14
// Not used
// double myapp::distance( const Eigen::Vector3d v3_pos1, const Eigen::Vector3d v3_pos2 ){
//   int i;
//   double sum2 = 0.0;
//   for(i=0;i<3;i++){
//     sum2 += (v3_pos2[i]-v3_pos1[i])*(v3_pos2[i]-v3_pos1[i]);
//   }
//   return sqrt(sum2);
// }


// 2022-06-06 11:22:37
// load index of nearest xy point
// of std::vector<Eigen::Vector3d> vec_v3;
// 2024-09-20 16:17:22 Not used.
// int myapp::load_index_of_nearest_xy( const std::vector<Eigen::Vector3d> &vec_XYZ, const Eigen::Vector3d &v3_pos ){
//   int i, index = -1;
//   double dx, dy, r2min = 9.9E+99;
//   for(i=0;i<vec_XYZ.size();i++){
//     // Eigen::Vector3d &xyz = vec_XYZ.at(i); // ERROR !
//     Eigen::Vector3d xyz = vec_XYZ.at(i);
//     dx = v3_pos.x()-xyz.x();
//     dy = v3_pos.y()-xyz.y();
//     if( r2min < dx*dx + dy*dy ){
//       r2min = dx*dx + dy*dy;
//       index = i;
//     }
//   }
//   assert(index != -1);
//   return index;
// }

// 2023-01-24 15:11:52
std::vector<std::pair<double,double>> myapp::read_double2( const std::filesystem::path &path_in )
{
  std::vector<std::pair<double,double>> vec_d_pair; // to be returned
  vec_d_pair.clear();  // Initialize

  std::ifstream reading_file;
  if(!std::filesystem::exists(path_in)) THROW_ERROR( path_in.string() + " does not exist");
  reading_file.open(path_in,std::ios::in);
  std::string reading_line_str;
  int nrow=0;
  while(std::getline(reading_file, reading_line_str)){
    // std::cerr << reading_line_str << std::endl;
    std::vector<std::string> vec_str = myapp::split(reading_line_str);
    if( vec_str.at(0).c_str()[0]==myapp::char_commentout_default ) continue;
    if( vec_str.size()!=2 ) THROW_ERROR("vec_str.size()!=2");
    double first  = std::stod( vec_str.at(0) );
    double second = std::stod( vec_str.at(1) );
    vec_d_pair.push_back( std::make_pair(first,second) );
    nrow++;
  }
  LOG_INFO("read {} rows", nrow);
  return vec_d_pair;
}

// for read 3clm ascii data
std::vector<Eigen::Vector3d> myapp::read_xyz( const std::filesystem::path &path_in )
{
  std::vector<Eigen::Vector3d> vec_v3; // Container
  vec_v3.clear();  // Initialize
  FILE *fin;
  char buf[1024];
  double x,y,z;
  Eigen::Vector3d v3(.0,.0,.0);
  int counter = 0;
  if(( fin=fopen(path_in.string().c_str(),"rt")) != NULL){
    while( fgets(buf,sizeof(buf),fin) != NULL ){
      if((buf[0] == myapp::char_commentout_default) || (buf[0] == '\n')) continue;
      if(counter%1000==0) fprintf(stderr,"counter = %d\r",counter);
      sscanf(buf,"%lf %lf %lf",&x,&y,&z);
      v3.x() = x;
      v3.y() = y;
      v3.z() = z;
      vec_v3.push_back(v3);
      counter++;
    }
  }
  else{ THROW_ERROR(" can't open file"); }
  if(fclose(fin) == EOF) THROW_ERROR("fclose(fin) == EOF, fname=" + path_in.string());
  LOG_INFO("# of point = {}.",vec_v3.size());
  
  return vec_v3;
}

// for read 3clm ascii data
std::vector<Eigen::Vector3d> myapp::read_xyz_stream( const std::filesystem::path &path_in )
{
  std::vector<Eigen::Vector3d> vec_v3; // to be returned
  vec_v3.clear();  // Initialize
  Eigen::Vector3d v3(.0,.0,.0);

  std::ifstream reading_file;
  if(!std::filesystem::exists(path_in)) THROW_ERROR( path_in.string() + " does not exist");
  reading_file.open(path_in,std::ios::in);
  std::string reading_line_str;
  int nrow=0;
  while(std::getline(reading_file, reading_line_str)){
    // std::cerr << reading_line_str << std::endl;
    // get divided vector string from line string
    std::vector<std::string> vec_str = myapp::split(reading_line_str);
    if( vec_str.at(0).c_str()[0]==myapp::char_commentout_default ) continue;
    if( vec_str.size()!=3 ) THROW_ERROR2("vec_str.size()!=3",vec_str.size());
    v3.x() = std::stod( vec_str.at(0) );
    v3.y() = std::stod( vec_str.at(1) );
    v3.z() = std::stod( vec_str.at(2) );
    vec_v3.push_back(v3);
    if(nrow%1000==0) fprintf(stderr,"counter = %d.....\r",nrow);
    nrow++;
  }
  LOG_INFO("read {} rows", nrow);
  return vec_v3;
}

// Read 3-column ASCII xyz data. Binary point-cloud input is no longer supported (throws).
std::vector< std::array<double,3> >
  myapp::read_vec_xyz( const std::filesystem::path &path_in )
{
  // path_in check
  myapp::filecheck(path_in);

  // Binary point-cloud input is no longer supported; only ASCII xyz is read here.
  // g2zbin DEM input is handled by the dedicated g2zbin reader, not this path.
  if (io_binary::is_binary_file(path_in)) {
    THROW_ERROR("myapp::read_vec_xyz: binary point-cloud input is no longer supported. path={}", path_in.string());
  }
  LOG_INFO("{} is text file", path_in.string());
  return myapp::read_vec_xyz_txt(path_in); // ascii file
}

// for read 3clm ascii data
std::vector< std::array<double,3> >
  myapp::read_vec_xyz_txt( const std::filesystem::path &path_in )
{
  std::vector<std::array<double,3>> vec_xyz; // to be returned

  std::ifstream reading_file;
  if(!std::filesystem::exists(path_in)) THROW_ERROR( path_in.string() + " does not exist");
  reading_file.open(path_in,std::ios::in);
  std::string reading_line_str;
  int nrow=0;
  double x,y,z;
  constexpr int EXPECTED_NUMBER_OF_VALUES = 3;
  constexpr int PRINT_INTERVAL = 1000;
  while(std::getline(reading_file, reading_line_str)){

    // Read entire line string, split by spaces, and put in vec_str.
    auto vec_str = myapp::split(reading_line_str);

    // Skip empty lines
    if (vec_str.empty()) {
      std::cerr << "Empty line at row=" + std::to_string(nrow+1) << std::endl;
      continue;
    }
    // Skip comment lines
    if( vec_str.at(0).c_str()[0]==myapp::char_commentout_default ) continue;

    // Error if not 3 columns of data
    if( vec_str.size()!=EXPECTED_NUMBER_OF_VALUES ) THROW_ERROR2("vec_str.size()!=3",vec_str.size());

    // Convert 3 columns of data from string to double
    x = std::stod( vec_str.at(0) );
    y = std::stod( vec_str.at(1) );
    z = std::stod( vec_str.at(2) );

    // Put 3 columns of data into vec_tp_xyz
    vec_xyz.push_back({x,y,z});

    if(nrow%PRINT_INTERVAL==0) fprintf(stderr,"counter = %d.....\r",nrow);
    nrow++;
  }
  fprintf(stderr,"\nread %d rows\n",nrow);
  if(nrow!=vec_xyz.size()){
    THROW_ERROR3("nrow!=vec_xyz.size()",nrow,vec_xyz.size());
  }
  return vec_xyz;
}

// for read 3double + 3 int ascii data
std::vector< std::tuple<double,double,double,int,int,int> >
  myapp::read_vec_xyzrgb( const std::filesystem::path &path_in )
{
  std::vector<std::tuple<double,double,double,int,int,int>> vec_tp_xyzrgb; // to be returned

  std::ifstream reading_file;
  if(!std::filesystem::exists(path_in)) THROW_ERROR( path_in.string() + " does not exist");
  reading_file.open(path_in,std::ios::in);
  std::string reading_line_str;
  int nrow=0;
  double x,y,z;
  int r,g,b;
  constexpr int EXPECTED_NUMBER_OF_VALUES = 6;
  constexpr int PRINT_INTERVAL = 1000;
  while(std::getline(reading_file, reading_line_str)){

    // Read entire line string, split by spaces, and put in vec_str.
    auto vec_str = myapp::split(reading_line_str);

    // Skip empty lines
    if (vec_str.empty()) {
      std::cerr << "Empty line at row=" + std::to_string(nrow+1) << std::endl;
      continue;
    }
    // Skip comment lines
    if( vec_str.at(0).c_str()[0]==myapp::char_commentout_default ) continue;

    // Error if not 3 columns of data
    if( vec_str.size()!=EXPECTED_NUMBER_OF_VALUES ) THROW_ERROR2("vec_str.size()!=EXPECTED_NUMBER_OF_VALUES",vec_str.size());

    // Convert 3 columns of data from string to double
    x = std::stod( vec_str.at(0) );
    y = std::stod( vec_str.at(1) );
    z = std::stod( vec_str.at(2) );
    r = std::stoi( vec_str.at(3) );
    g = std::stoi( vec_str.at(4) );
    b = std::stoi( vec_str.at(5) );

    // Put 3 columns of data into vec_tp_xyz
    vec_tp_xyzrgb.push_back(std::make_tuple(x,y,z,r,g,b));

    if(nrow%PRINT_INTERVAL==0) fprintf(stderr,"counter = %d.....\r",nrow);
    nrow++;
  }
  fprintf(stderr,"\nread %d rows\n",nrow);
  if(nrow!=vec_tp_xyzrgb.size()){
    THROW_ERROR3("nrow!=vec_tp_xyz.size()",nrow,vec_tp_xyzrgb.size());
  }
  return vec_tp_xyzrgb;
}

bool myapp::compare_tuple_yxz(
  const std::tuple<double, double, double>& a,
  const std::tuple<double, double, double>& b)
{
  // Compare second element
  if (std::get<1>(a) < std::get<1>(b)) {
    return true;
  }
  else if (std::get<1>(a) > std::get<1>(b)) {
    return false;
  }
  else {
    // If second element is same, compare first element
    if (std::get<0>(a) < std::get<0>(b)) {
      return true;
    }
    else if (std::get<0>(a) > std::get<0>(b)) {
      return false;
    }
    else {
      // If second and first elements are same, compare third element
      if (std::get<2>(a) < std::get<2>(b)) {
        return true;
      }
      else {
        return false;
      }
    }
  }
}

// @brief for sort. compare array<double,3> by yxz
bool myapp::compare_array_yxz(
  const std::array<double, 3>& a, const std::array<double, 3>& b)
{
  // Compare second element
  if (a[1] < b[1]) {
    return true;
  }
  else if (a[1] > b[1]) {
    return false;
  }
  else {
    // If second element is same, compare first element
    if (a[0] < b[0]) {
      return true;
    }
    else if (a[0] > b[0]) {
      return false;
    }
    else {
      // If second and first elements are same, compare third element
      if (a[2] < b[2]) {
        return true;
      }
      else {
        return false;
      }
    }
  }
}

// output std::vector<Eigen::Vector3d>
void myapp::out_vec_v3( const std::filesystem::path& pathout, const std::vector<Eigen::Vector3d> &vec_v3 )
{
  FILE *fout = get_fout(pathout);
  // for(const auto& v : vec_v3) fprintf(fout,"%E %E %E\n",v.x(),v.y(),v.z());
  for(const auto& v : vec_v3) fprintf(fout,"%E %E %E\n",v.x(),v.y(),v.z());
  myapp::close(fout,pathout);
}

//   std::vector< std::array<double,3> >
void myapp::out( const std::filesystem::path& pathout
  , const std::vector< std::array<double,3>> &vec_tp )
{
  FILE *fout = get_fout(pathout);
  // for(auto v : vec_v3) fprintf(fout,"%E %E %E\n",v.x(),v.y(),v.z());
  for(const auto& [x,y,z] : vec_tp) fprintf(fout,"%E %E %E\n",x,y,z);
  myapp::close(fout,pathout);
}

//   std::vector< std::array<double,4> >
void myapp::out( const std::filesystem::path& pathout
  , const std::vector< std::array<double,4>> &vec_tp4 )
{
  FILE *fout = get_fout(pathout);
  // for(auto v : vec_v3) fprintf(fout,"%E %E %E\n",v.x(),v.y(),v.z());
  for(const auto& [x,y,z,w] : vec_tp4) fprintf(fout,"%E %E %E %E\n",x,y,z,w);
  myapp::close(fout,pathout);
}

//   std::vector< std::array<double,4> >
void myapp::out_vec_xyzdzrgb( const std::filesystem::path& pathout
  , const std::vector< std::tuple<double,double,double,double,int,int,int>> &vec_xyzdzrgb )
{
  FILE *fout = get_fout(pathout);
  // for(auto v : vec_v3) fprintf(fout,"%E %E %E\n",v.x(),v.y(),v.z());
  for(const auto& [x,y,z,w,r,g,b] : vec_xyzdzrgb) fprintf(fout,"%E %E %E %E %d %d %d\n",x,y,z,w,r,g,b);
  myapp::close(fout,pathout);
}

// 2022-09-12 10:35:11
// print all std::vector<std::pair<int,int>>
void myapp::out_vec_pair( FILE *fout, const std::vector<std::pair<int,int>> &vec_ipair )
{
  for( const auto& v : vec_ipair ) fprintf(fout,"%5d %5d\n",v.first,v.second);
}

void myapp::out_vec_pair(
  const std::filesystem::path& pathout, const std::vector<std::pair<int,int>> &vec_ipair )
{
  FILE *fout = get_fout(pathout);
  out_vec_pair(fout,vec_ipair);
  myapp::close(fout,pathout);
}

// 2022-09-16 13:33:01
// print all std::vector<std::pair<int,int>>
void myapp::out_vec_pair( FILE *fout, const std::vector<std::pair<double,double>> &vec_pair )
{
  for(const auto& [x, y] : vec_pair) fprintf(fout, "%E %E\n", x, y);
}

void myapp::out_vec_pair(
  const std::filesystem::path& pathout, const std::vector<std::pair<double,double>> &vec_pair )
{
  FILE *fout = get_fout(pathout);
  out_vec_pair(fout,vec_pair);
  myapp::close(fout,pathout);
}

// 2022-09-16 14:51:34
// for debug
void myapp::debug_y_ax_b(FILE *fout, const std::vector<std::pair<double,double>> &vec_pair, const double &a, const double &b ){
  for(const auto& v : vec_pair ){
    double y = a * v.first + b;
    double dy = v.second - y;
    fprintf(fout,"%E %E %E %E\n",v.first,v.second,y,dy);
  }
}

void myapp::debug_y_ax_b(const std::filesystem::path& pathout
, const std::vector<std::pair<double,double>> &vec_pair, const double &a, const double &b )
{
  FILE *fout = get_fout(pathout);
  debug_y_ax_b(fout,vec_pair,a,b);
  myapp::close(fout,pathout);
}

// for output std::vector<std::tuple<int,int,int>> to file
void myapp::out_int_triple(FILE *fout
, const std::tuple<int,int,int> &tp)
{
  const int& ix = std::get<0>(tp);
  const int& iy = std::get<1>(tp);
  const int& iz = std::get<2>(tp);
  fprintf(fout,"%5d %5d %5d\n",ix,iy,iz);
}

void myapp::out_vec_triple( FILE *fout
, const std::vector<std::tuple<int,int,int>> &vec_triple )
{
  for(const auto& tp : vec_triple ) out_int_triple(fout,tp);
}

void myapp::out_vec_triple( const std::filesystem::path& pathout
, const std::vector<std::tuple<int,int,int>> &vec_triple )
{
  FILE *fout = get_fout(pathout);
  out_vec_triple(fout,vec_triple);
  myapp::close(fout,pathout);
}

// for output std::vector<std::tuple<int,int,double>> to file
void myapp::out_ix_iy_PL(FILE *fout
, const std::tuple<int,int,double> &tp)
{
  const int& ix = std::get<0>(tp);
  const int& iy = std::get<1>(tp);
  const double& path = std::get<2>(tp);
  fprintf(fout,"%5d %5d %lf\n",ix,iy,path);
}

void myapp::out_vec_ix_iy_PL( FILE *fout
, const std::vector<std::tuple<int,int,double>> &vec_ix_iy_PL )
{
  for( const auto& tp : vec_ix_iy_PL ) out_ix_iy_PL(fout,tp);
}

void myapp::out_vec_ix_iy_PL( const std::filesystem::path& pathout
, const std::vector<std::tuple<int,int,double>> &vec_ix_iy_PL )
{
  FILE *fout = get_fout(pathout);
  out_vec_ix_iy_PL(fout,vec_ix_iy_PL);
  myapp::close(fout,pathout);
}

// for output std::vector<std::tuple<int,int,double>> to file
void myapp::out_ix_iy_iz_path(FILE *fout
, const std::tuple<int,int,int,double> &tp)
{
  const int& ix = std::get<0>(tp);
  const int& iy = std::get<1>(tp);
  const int& iz = std::get<2>(tp);
  const double& path = std::get<3>(tp);
  fprintf(fout,"%5d %5d %5d %lf\n",ix,iy,iz,path);
}

void myapp::out_vec_ix_iy_iz_path( FILE *fout
, const std::vector<std::tuple<int,int,int,double>> &vec_ix_iy_iz_path )
{
  for(const auto& tp : vec_ix_iy_iz_path ) out_ix_iy_iz_path(fout,tp);
}

void myapp::out_vec_ix_iy_iz_path( const std::filesystem::path& pathout
, const std::vector<std::tuple<int,int,int,double>> &vec_ix_iy_iz_path )
{
  FILE *fout = get_fout(pathout);
  out_vec_ix_iy_iz_path(fout,vec_ix_iy_iz_path);
  myapp::close(fout,pathout);
}

/// @brief erase all files in directory. directory is not erased.
void myapp::clear_dir(const std::filesystem::path& dir_path)
{
  LOG_INFO("dir_path = {}", dir_path.string());
  LOG_INFO("this function erases all files in the directory, but not the directory itself.");
  LOG_INFO("this function also erases subdirectories and their contents.");
  namespace fs = std::filesystem;
  if (!fs::exists(dir_path) || !fs::is_directory(dir_path)) {
    LOG_ERROR("Invalid path: either the directory does not exist or the path is not a directory.");
    LOG_WARN("please check also the order to apply this function.");
    LOG_WARN("e.g. std::filesystem::current_path({}); then myapp::clear({}) generates error.");
    LOG_WARN("In this case, currect directory {} should be above the directory you want to clear.", dir_path.string());
    THROW_ERROR("myapp::clear_dir Invalid path: either the directory does not exist or the path is not a directory.");
  }

  for (const auto& entry : fs::directory_iterator(dir_path)) {
    try {
      fs::remove_all(entry.path());  // OK whether it is a file or a subdirectory
    } catch (const std::exception& e) {
      LOG_ERROR("Deletion failed: {}: {}", entry.path().string(), e.what());
    }
  }
}

// time measurement
// 2022-09-30 17:21:54
void myapp::cast_time_min( FILE *fout
  , const std::string &msg
  , const std::chrono::system_clock::time_point &start
  , const std::chrono::system_clock::time_point &end )
{
  const double time_sec = 
    static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(end - start).count());
  fprintf(fout,"##################################################################\n");
  fprintf(fout,"time = %.0lf min, %s\n", time_sec/60,msg.c_str());
  fprintf(fout,"##################################################################\n");
}

void myapp::cast_time_min(
    const spdlog::level::level_enum log_level
  , const std::string &msg
  , const std::chrono::system_clock::time_point &start
  , const std::chrono::system_clock::time_point &end )
{
  const double time_sec = static_cast<double>(
    std::chrono::duration_cast<std::chrono::seconds>(end - start).count());

  const std::string decorated_msg = fmt::format(
    "\n"
    "##################################################################\n"
    "time = {:.0f} min, {}\n"
    "##################################################################"
  , time_sec / 60.0, msg);

  mylogger::safe_log(log_level, decorated_msg);
}

void myapp::cast_time_sec( FILE *fout
  , const std::string &msg
  , const std::chrono::system_clock::time_point &start
  , const std::chrono::system_clock::time_point &end )
{
  const double time_sec = 
    static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(end - start).count());
  fprintf(fout,"##################################################################\n");
  fprintf(fout,"time =  %.0lf sec, %s\n",time_sec,msg.c_str());
  fprintf(fout,"##################################################################\n");
}

void myapp::cast_time_sec(
    const spdlog::level::level_enum log_level
  , const std::string &msg
  , const std::chrono::system_clock::time_point &start
  , const std::chrono::system_clock::time_point &end )
{
  const double time_sec = static_cast<double>(
    std::chrono::duration_cast<std::chrono::seconds>(end - start).count());

  const std::string decorated_msg = fmt::format(
    "\n"
    "##################################################################\n"
    "time = {:.0f} sec, {}\n"
    "##################################################################"
  , time_sec, msg);

  mylogger::safe_log(log_level, decorated_msg);
}

void myapp::cast_time_msec( FILE *fout
  , const std::string &msg
  , const std::chrono::system_clock::time_point &start
  , const std::chrono::system_clock::time_point &end )
{
  const double time_millisec = 
    static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
  fprintf(fout,"##################################################################\n");
  fprintf(fout,"time = %.0lf msec, %s\n",time_millisec, msg.c_str());
  fprintf(fout,"##################################################################\n");
}

void myapp::cast_time_msec(
    const spdlog::level::level_enum log_level
  , const std::string &msg
  , const std::chrono::system_clock::time_point &start
  , const std::chrono::system_clock::time_point &end )
{
  const double time_millisec = static_cast<double>(
    std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

  const std::string decorated_msg = fmt::format(
    "\n"
    "##################################################################\n"
    "time = {:.0f} msec, {}\n"
    "##################################################################"
  , time_millisec, msg);

  mylogger::safe_log(log_level, decorated_msg);
}



// 2022-11-25 11:57:24
// out_vec_double
void myapp::out( const std::filesystem::path& pathout, const std::vector<double> &vec )
{
  FILE *fout = get_fout(pathout);
  for(const auto& val : vec ) fprintf(fout,"%E\n",val);
  myapp::close(fout,pathout);
  LOG_INFO("std::vector<double> {} done.",pathout.string());
}


//2023-03-29 12:58:29
void myapp::out( const std::filesystem::path& pathout, const Eigen::VectorXf &vecxf )
{
  FILE *fout = get_fout(pathout);
  const int nrow = vecxf.rows();
  for(int irow=0;irow<nrow;irow++) fprintf(fout,"%d %E\n",irow,vecxf(irow));
  myapp::close(fout,pathout);
  LOG_INFO("Eigen::VectorXf {} done.",pathout.string());
}

// mabe by chatGPT 3.5
// check interval x
bool myapp::is_same_interval_x(
  const std::vector<std::pair<double,double>> &vec_double2, const double double_epsilon )
{
  if (vec_double2.size() < 2) {
    return false; // An empty vector or a single-element vector is not considered evenly spaced
  }
  const double diff = vec_double2.at(1).first - vec_double2.at(0).first; // First difference
  for (size_t i = 2; i < vec_double2.size(); i++) {
    const double cur_diff = vec_double2.at(i).first - vec_double2.at(i-1).first;
    const double abs_diff = std::abs( std::abs(cur_diff) - std::abs(diff) );
    const double abs_sum = std::abs( std::abs(cur_diff) + std::abs(diff) ) * double_epsilon;
    if ( abs_diff > abs_sum ) {
      return false; // Not evenly spaced if the difference is not constant
    }
  }
  return true; // Evenly spaced if the difference is constant
}

//===================================================
// made my chatGPT
// enum class VectorOrder
// 2023-04-07 13:45:53
//===================================================
VectorOrder myapp::getVectorOrder(const std::vector<double>& vec)
{
  bool isIncreasing = true;
  bool isDecreasing = true;

  for (int i = 1; i < vec.size(); ++i) {
    if (vec.at(i) > vec.at(i-1)) {
      isDecreasing = false;
    } else if (vec.at(i) < vec.at(i-1)) {
      isIncreasing = false;
    }
  }

  if (isIncreasing) {
    return VectorOrder::Increasing;
  } else if (isDecreasing) {
    return VectorOrder::Decreasing;
  } else {
    return VectorOrder::NotOrdered;
  }
}

// 2023-04-14 14:41:43
// by chatGPT
MonotonicStatus myapp::getMonotonicStatus(const std::vector<std::pair<double, double>>& data)
{
  if (data.size() < 2) {
    return MonotonicStatus::NotMonotonic;
  }

  int direction = 0;  // 0: undetermined, 1: increasing, -1: decreasing

  for (size_t i = 1; i < data.size(); ++i) {
    double dx = data[i].first - data[i - 1].first;
    double dy = data[i].second - data[i - 1].second;

    if (dx == 0 || dy == 0) {
      continue;
    }

    int current_direction;
    if (dx > 0 && dy > 0) {
      current_direction = 1;
    } else if (dx < 0 && dy < 0) {
      current_direction = 1;
    } else {
      current_direction = -1;
    }
    if (direction == 0) {
      direction = current_direction;
    } else if (direction != current_direction) {
      return MonotonicStatus::NotMonotonic;
    }
  }

  MonotonicStatus status;
  if (direction == 0) {
    status = MonotonicStatus::NotMonotonic;
  } else if (direction == 1) {
    status = MonotonicStatus::Increasing;
  } else {
    status = MonotonicStatus::Decreasing;
  }
  return status;
}


//====================================================
// Interpolation calculation functions
//====================================================
double myapp::linear_interpolation(double x0, double y0, double x1, double y1, double x) {
  // Tolerance to absorb floating-point rounding when x lands exactly on a node,
  // matching the guard already used by bilinear_interpolation below. Without it,
  // a query sitting one ULP past the upper node (a node-aligned, in-range value)
  // is wrongly rejected; with it, the formula returns the node value as expected.
  constexpr double small = 1000.0 * std::numeric_limits<double>::epsilon();
  if (x0 < x1) {
    if (x + small < x0) THROW_ERROR("angle_util::linear_interpolation x  < x0");
    if (x1 + small <= x) THROW_ERROR("angle_util::linear_interpolation x1 <= x");
  } else if (x1 < x0) {
    if (x + small < x1) THROW_ERROR("angle_util::linear_interpolation x  <  x1");
    if (x0 + small <= x) THROW_ERROR("angle_util::linear_interpolation x0 <= x ");
  } else {
    THROW_ERROR("angle_util::linear_interpolation x0==x1");
  }
  return (y1 - y0) / (x1 - x0) * (x - x0) + y0;
}

double myapp::bilinear_interpolation(double x1, double x2, double y1, double y2,
                              double f11, double f21, double f12, double f22,
                              double x, double y) {
  constexpr double small = 1000.0 * std::numeric_limits<double>::epsilon();
  if (x2 + small < x1) THROW_ERROR("angle_util::bilinear_interpolation,x2 + small <  x1");
  if (x + small < x1) THROW_ERROR("angle_util::bilinear_interpolation,x  + small <  x1");
  if (x2 + small <= x) THROW_ERROR("angle_util::bilinear_interpolation,x2 + small <= x ");
  if (y2 + small < y1) THROW_ERROR("angle_util::bilinear_interpolation,y2 + small <  y1");
  if (y + small < y1) THROW_ERROR("angle_util::bilinear_interpolation,y  + small <  y1");
  if (y2 + small <= y) THROW_ERROR("angle_util::bilinear_interpolation,y2 + small <= y ");
  const double norm_factor = 1 / ((x2 - x1) * (y2 - y1));
  const double term11 = f11 * (x2 - x) * (y2 - y);
  const double term21 = f21 * (x - x1) * (y2 - y);
  const double term12 = f12 * (x2 - x) * (y - y1);
  const double term22 = f22 * (x - x1) * (y - y1);
  return norm_factor * (term11 + term21 + term12 + term22);
}

//==================================================
// std::vector< std::array<double,3> > 
// Function that checks whether they match
//==================================================
// Comparison of primitive floating-point numbers (allowed tolerance specified via epsilon)
bool myapp::is_eq( const double a, const double b, const double epsilon)
{
  if( a==0.0 && b==0.0 ) return true;
  const double small = epsilon*std::abs(a + b);
  const double diff = std::abs(a - b);
  return diff < small;
}

// Comparison of std::tuple<double, double, double>
bool myapp::is_eq(
    const std::tuple<double, double, double>& t1
  , const std::tuple<double, double, double>& t2
  , const double epsilon)
{
  const auto [x1,y1,z1] = t1;
  const auto [x2,y2,z2] = t2;
  bool tf_x = myapp::is_eq(x1,x2,epsilon);
  if( !tf_x ) return false;
  
  bool tf_y = myapp::is_eq(y1,y2,epsilon);
  if( !tf_y ) return false;
  
  bool tf_z = myapp::is_eq(z1,z2,epsilon);
  if( !tf_z ) return false;

  return true;
}

// Comparison of std::vector< std::tuple<double, double, double> >
bool myapp::is_eq(
    const std::vector< std::tuple<double, double, double> >& v1
  , const std::vector< std::tuple<double, double, double> >& v2
  , const double epsilon)
{
  if (v1.size() != v2.size()) return false;

  for (size_t i = 0; i < v1.size(); ++i) {
    const auto tp1 = v1.at(i);
    const auto tp2 = v2.at(i);
    if (!myapp::is_eq(tp1, tp2, epsilon)) {
      return false;
    }
  }
  return true;
}

// Splits the arr array into N groups.
// Assigns the remainder to the targetGroup-th group.
// If targetGroup=0, the remainder is not assigned.
std::vector<std::vector<int>> myapp::splitArrayIntoN(
  const std::vector<int>& arr, int N, const int targetGroup)
{
  if( N>=arr.size() ) THROW_ERROR("N>=arr.size()");
  if( targetGroup>N ) THROW_ERROR("targetGroup>N");
  std::vector<std::vector<int>> groups(N);
  int size = arr.size();
  int baseLength = size / N;
  int extra = size % N;

  int start = 0;

  for (int i = 0; i < N; ++i) {
    int length = baseLength;

    if (i == targetGroup -1 ) {
      length += extra;
    }
    
    for (int j = start; j < start + length; ++j) {
      groups.at(i).push_back(arr.at(j));
    }
    start += length;
  }

  return groups;
}

// Splits the arr array into nx_div x ny_div groups.
// The remainder is assigned to targetGroup_x and targetGroup_y.
std::vector<std::vector<std::vector<std::vector<std::array<int,2>>>>> 
  myapp::split2DArrayIntoGrid(
    const std::vector<std::vector<std::array<int,2>>>& arr
  , const int nx_div, const int ny_div, const int targetGroup_x, const int targetGroup_y)
{
  int ny = arr.size();
  int nx = arr.at(0).size();

  std::vector<std::vector<std::vector<std::vector<std::array<int,2>>>>> result(ny_div, std::vector<std::vector<std::vector<std::array<int,2>>>>(nx_div));

  int baseRowLength = ny / ny_div;
  int extraRow = ny % ny_div;

  int baseColLength = nx / nx_div;
  int extraCol = nx % nx_div;

  int startRow = 0;
  for (int iy = 0; iy < ny_div; ++iy) {
    int rowLength = baseRowLength;
    if (iy == targetGroup_y-1) {
      rowLength += extraRow;
    }

    int startCol = 0;
    for (int ix = 0; ix < nx_div; ++ix) {
      int colLength = baseColLength;
      if (ix == targetGroup_x-1) {
        colLength += extraCol;
      }

      for (int y = startRow; y < startRow + rowLength; ++y) {
        std::vector<std::array<int,2>> temp;
        for (int x = startCol; x < startCol + colLength; ++x) {
          temp.push_back(arr.at(y).at(x));
        }
        result[iy][ix].push_back(temp);
      }

      startCol += colLength;
    }

    startRow += rowLength;
  }

  return result;
}

void myapp::printResult(
  const std::vector<std::vector<std::vector<std::vector<std::array<int,2>>>>>& result)
{
  // Outer loop is vertical direction of Group
  for (int i = 0; i < result.size(); ++i) {

    // Inner loop is horizontal direction of Group
    for (int j = 0; j < result[i].size(); ++j) {
      
      // Display summary
      std::cout << "Group(" << i + 1 << "," << j + 1 << ") has "
                << result[i][j].size() << " rows "
                << (result[i][j].empty() ? 0 : result[i][j][0].size()) << " columns." << std::endl;

      // Display Group body
      for (const auto& row : result[i][j]) {
        for (const auto& elem : row) {
          const auto [ix, iy] = elem;
          printf(" (%02d,%02d)", ix, iy);
        }
        std::cout << std::endl;
      }
      std::cout << std::endl; // Insert blank line between Groups
    }
  }
}

void myapp::printGroups(const std::vector<std::vector<std::vector<std::array<int,2>>>>& groups)
{
  int ny_div = groups.size();
  int nx_div = groups[0].size();

  for (int iy = 0; iy < ny_div; ++iy) {
    for (int ix = 0; ix < nx_div; ++ix) {
      int row = groups[iy][ix].size();
      std::cout << "Group(" << iy + 1 << "," << ix + 1 << ") has " << row << " rows.\n";
      
      for (const auto& r : groups[iy][ix]) {
        const auto [i,j] = r;
        std::cout << "(" << i << "," << j << ") ";
      }
      std::cout << "\n";
    }
  }
}

// Function to copy int corresponding to specified num from std::vector<std::array<int,3>>
// to std::vector<int>
std::vector<int> myapp::copyArrayElementToVector(
  const std::vector<std::array<int,3>>& input, const int num )
{
  std::vector<int> output;

  for (const auto& elem : input) {
    if (num >= 0 && num < 3) {
      output.push_back(elem[num]);
    } else {
      THROW_ERROR("myapp::copyArrayElementToVector: Invalid num value. num={}, valid range=[0,2]", num);
    }
  }

  return output;
}

// In this version, std::string object is initialized with a pre-specified size.
// Then, the region is directly read using ifs.read method.
// This method is efficient, but can be said to depend on the internal implementation of std::string
// (though this is legal from C++11 onwards).


// Helper function for reading
//
// std::string read_string(std::ifstream& ifs) {
//   int length;
//   ifs.read(reinterpret_cast<char*>(&length), sizeof(int));
//   char* buffer = new char[length];
//   ifs.read(buffer, length);
//   std::string str(buffer, length);
//   delete[] buffer;
//   return str;
// }
// In this version, reading is done into a temporary buffer and then copied to std::string.
// This method is slightly redundant, but has the advantage of not depending on the internal representation of std::string.






// Function to determine whether two tuple_int::UmpIntInt3 are equal
bool myapp::are_ump_eq(const tuple_int::UmpIntInt3& map1, const tuple_int::UmpIntInt3& map2 )
{
  // If sizes differ, maps do not match
  if (map1.size() != map2.size()) {
    return false;
  }

  // Check if each element of map1 also exists in map2
  for (const auto& [key, value] : map1) {
    auto it = map2.find(key);
    if (it == map2.end()) {
      return false; // Key does not exist
    }

    if (it->second != value) {
      return false; // Value does not match
    }
  }

  return true;
}

// Function to determine whether two tuple_int::UmpInt3Int are equal
bool myapp::are_ump_eq(const tuple_int::UmpInt3Int& map1, const tuple_int::UmpInt3Int& map2 )
{
  // If sizes differ, maps do not match
  if (map1.size() != map2.size()) {
    return false;
  }

  // Check if each element of map1 also exists in map2
  for (const auto& [key, value] : map1) {
    auto it = map2.find(key);
    if (it == map2.end()) {
      return false; // Key does not exist
    }

    if (it->second != value) {
      return false; // Value does not match
    }
  }

  return true;
}


// Function to check continuity of vector<int>.
// Empty vector is considered continuous.
bool myapp::sort_and_check_vec_int_continuity( std::vector<int>& vec_int )
{
  if (vec_int.empty()) return true; // Empty vector is considered continuous

  std::sort(vec_int.begin(), vec_int.end());

  const int min = vec_int.front();
  const int max = vec_int.back();
  return ((max - min + 1) == vec_int.size());
}

void myapp::out_set( const std::set<int>& data, const std::filesystem::path& pathout)
{
  std::ofstream ofs(pathout);
  if (!ofs.is_open()) {
    LOG_ERROR("write_set_to_file: failed to open file {}", pathout.string());
    THROW_ERROR(fmt::format("write_set_to_file: cannot open file {}", pathout.string()));
  }

  for (const auto& val : data) {
    ofs << val << "\n";
  }

  LOG_INFO("write_set_to_file: wrote {} entries to {}", data.size(), pathout.string());
}

// for output std::map<int, bool>
void myapp::out_map(const std::map<int, bool>& map_int_bool, const std::filesystem::path& pathout)
{
  std::ofstream file(pathout);
  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << pathout <<  std::endl;
    return;
  }

  for (const auto& pair : map_int_bool) {
    file << pair.first << " " << (pair.second ? "true" : "false") << std::endl;
  }

  file.close();
}

void myapp::out_map(const std::map<int, int>& map_int_int, const std::filesystem::path& pathout)
{
  std::ofstream file(pathout);
  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << pathout <<  std::endl;
    return;
  }
  for (const auto& pair : map_int_int) {
    file << pair.first << " " << pair.second << std::endl;
  }

  file.close();
}

// for output std::unordered_multimap<int, int>
void myapp::out_uommp_int2(const std::unordered_multimap<int,int>& uommp, const std::filesystem::path& pathout)
{
  std::ofstream outFile(pathout);
  if (!outFile.is_open()) {
    std::cerr << "Cannot open file: " << pathout << std::endl;
    return;
  }

  for (const auto& pair : uommp) {
    outFile << pair.first << " -> " << pair.second << '\n';
  }
  fprintf(stderr, "myapp::out_uommp_int2, File output: %s\n", pathout.string().c_str());
  outFile.close();
}

// for output std::unordered_map<int, int>
void myapp::out_uomp_int2(const std::unordered_map<int,int>& uomp
  , const std::filesystem::path& pathout)
{
  std::ofstream outFile(pathout);
  if (!outFile.is_open()) {
    std::cerr << "Cannot open file: " << pathout << std::endl;
    return;
  }

  for (const auto& pair : uomp) {
    outFile << pair.first << " -> " << pair.second << '\n';
  }
  fprintf(stderr, "myapp::out_uomp_int2, File output: %s\n", pathout.string().c_str());
  outFile.close();
}

void myapp::write_sorted_detid_uqigavail_to_file(
  const SortedDetidUqigSet& data
, const std::filesystem::path& pathout)
{
  std::ofstream ofs(pathout);
  if (!ofs.is_open()) {
    LOG_ERROR("write_sorted_detid_uqigavail_to_file: failed to open file {}", pathout.string());
    THROW_ERROR(fmt::format("write_sorted_detid_uqigavail_to_file: cannot open file {}", pathout.string()));
  }

  for (const auto& [detid, uqigavail] : data) {
    ofs << detid << "->" << uqigavail << "\n";
  }

  LOG_INFO("write_sorted_detid_uqigavail_to_file: wrote {} entries to {}"
  , data.size(), pathout.string());
}


// compare unordered_multimap<int,int> is same or not
bool myapp::is_eq(const std::unordered_multimap<int, int>& ummap1, const std::unordered_multimap<int, int>& ummap2)
{
  // Compare sizes
  if (ummap1.size() != ummap2.size()) {
    return false;
  }

  for (const auto& pair : ummap1) {
    auto range = ummap2.equal_range(pair.first);
    int count = ummap1.count(pair.first);

    if (std::distance(range.first, range.second) != count) {
      return false;
    }

    for (auto it = range.first; it != range.second; ++it) {
      if (it->second != pair.second) {
        return false;
      }
    }
  }

  return true;
}

//#####################################################################################
/// @namespace my_json_func
/// @brief Functions related to JSON file reading and formatted output
//#####################################################################################

// Function to load json file ignoring comments
nlohmann::json myapp::load_json(const std::filesystem::path &pathin)
{
  fprintf(stderr, "myapp::load_json, pathin = %s\n", pathin.string().c_str());
  std::ifstream file(pathin);
  if (!file.is_open()) {
    THROW_ERROR("myapp::load_json: Unable to open file. path={}", pathin.string());
  }
  std::stringstream buffer;
  buffer << file.rdbuf();
  std::string fileContents = buffer.str();
  // Pass true as 4th argument to ignore comments
  return nlohmann::json::parse(fileContents, nullptr, true, true);
}

// Function to format and output JSON
// usage1 : dump_with_section(js, {}, "address", spdlog::level::debug); --> output hierarchy up to "address" and "address" section
// usage2 : dump_with_section(js, {"profile", "address"}, "location", spdlog::level::debug); --> output only "location" section
void myapp::dump_json(
    const nlohmann::json& js
  , const std::vector<std::string>& current_path
  , const std::string& target_section
  , const spdlog::level::level_enum& log_level
  // Do not specify arguments below normally
  , const int indent
  , std::ostringstream* oss_ptr
  , bool already_wrapped
  , bool inline_mode /* = false */)
{
  // Wrap only at first (during recursion already_wrapped == true)
  nlohmann::json js_to_dump = (already_wrapped || js.is_array()) ? js : nlohmann::json::array({js});
  std::ostringstream oss_local;
  std::ostringstream& oss = oss_ptr ? *oss_ptr : oss_local;

  bool in_target = (target_section.empty() ||
    std::find(current_path.begin(), current_path.end(), target_section) != current_path.end());

  if (!in_target && !current_path.empty()) {
    oss << std::string(indent * current_path.size(), ' ') << "\"...\"\n";
    if (!oss_ptr)
      mylogger::g_logger->log(log_level, oss.str());
    return;
  }

  if (js_to_dump.is_object()) {
    oss << std::string(indent * current_path.size(), ' ') << "{\n";
    for (auto it = js_to_dump.begin(); it != js_to_dump.end(); ++it) {
      oss << std::string(indent * (current_path.size() + 1), ' ')
          << "\"" << it.key() << "\": ";
      auto next_path = current_path;
      next_path.push_back(it.key());
      // Inside object, set inline_mode = false to insert line breaks
      dump_json(it.value(), next_path, target_section, log_level, indent, &oss, true, false);
    }
    oss << std::string(indent * current_path.size(), ' ') << "}";
    if (!inline_mode)
      oss << "\n";
  } else if (js_to_dump.is_array()) {
    oss << std::string(indent * current_path.size(), ' ') << "[\n";
    size_t index = 0;
    const size_t size = js_to_dump.size();
    for (const auto& el : js_to_dump) {
      oss << std::string(indent * (current_path.size() + 1), ' ');
      // If array element is primitive, call with inline_mode = true
      if (!el.is_object() && !el.is_array()) {
        dump_json(el, current_path, target_section, log_level, indent, &oss, true, true);
      }
      else {
        dump_json(el, current_path, target_section, log_level, indent, &oss, true, false);
      }
      if (index < size - 1) {
        oss << ",";
      }
      oss << "\n";
      ++index;
    }
    oss << std::string(indent * current_path.size(), ' ') << "]";
    if (!inline_mode)
      oss << "\n";
  } else {
    // For primitive, output with dump() (no line break if inline_mode true)
    oss << js_to_dump.dump();
    if (!inline_mode) {
      oss << "\n";
    }
  }

  if (!oss_ptr)
    mylogger::g_logger->log(log_level, oss.str());
}

void myapp::dump_json(
    const nlohmann::json& js
  , const std::vector<std::string>& current_path
  , const std::string& target_section
  , FILE* fout
  // Do not specify arguments below normally
  , const int indent
  , bool already_wrapped
  , bool inline_mode /* = false */)
{
  // Wrap only at first (during recursion already_wrapped==true)
  nlohmann::json js_to_dump = (already_wrapped || js.is_array()) ? js : nlohmann::json::array({js});

  bool in_target = (target_section.empty() ||
    std::find(current_path.begin(), current_path.end(), target_section) != current_path.end());

  if (!in_target && !current_path.empty()) {
    fprintf(fout, "%*s\"...\"\n", static_cast<int>(indent * current_path.size()), "");
    return;
  }

  if (js_to_dump.is_object()) {
    fprintf(fout, "%*s{\n", static_cast<int>(indent * current_path.size()), "");
    size_t index = 0;
    const size_t size = js_to_dump.size();
    for (auto it = js_to_dump.begin(); it != js_to_dump.end(); ++it, ++index) {
      fprintf(fout, "%*s\"%s\": ", static_cast<int>(indent * (current_path.size() + 1)), "",
              it.key().c_str());
      auto next_path = current_path;
      next_path.push_back(it.key());
      // Inside object, insert line break (inline_mode=false)
      dump_json(it.value(), next_path, target_section, fout, indent, true, false);
      if (index < size - 1) {
        fprintf(fout, ",");
      }
      // fprintf(fout, "\n");
    }
    fprintf(fout, "%*s}", static_cast<int>(indent * current_path.size()), "");
    if (!inline_mode) fprintf(fout, "\n");
  }
  else if (js_to_dump.is_array()) {
    fprintf(fout, "%*s[\n", static_cast<int>(indent * current_path.size()), "");
    size_t index = 0;
    const size_t size = js_to_dump.size();
    for (const auto& el : js_to_dump) {
      fprintf(fout, "%*s", static_cast<int>(indent * (current_path.size() + 1)), "");
      // For primitive, call with inline_mode=true (same line output)
      if (!el.is_object() && !el.is_array()) {
        dump_json(el, current_path, target_section, fout, indent, true, true);
      }
      else {
        dump_json(el, current_path, target_section, fout, indent, true, false);
      }
      if (index < size - 1) {
        fprintf(fout, ",");
      }
      if (!inline_mode) fprintf(fout, "\n");
      ++index;
    }
    fprintf(fout, "%*s]", static_cast<int>(indent * current_path.size()), "");
    if (!inline_mode) fprintf(fout, "\n");
  }
  else {
    // For primitive, dump() output (no line break at end if inline_mode is true)
    fprintf(fout, "%s", js_to_dump.dump().c_str());
    if (!inline_mode) fprintf(fout, "\n");
  }
}

//#####################################################################################
/// @namespace mysort
//#####################################################################################

// Sort std::vector<std::pair<double,double>> in ascending order by first
void mysort::sort_vec_double2_x_increasing( std::vector<std::pair<double,double>> &vec_double2 )
{
  std::sort(vec_double2.begin(),vec_double2.end(),mysort::vec_double2_sort_x_increasing());
};

// Sort std::vector<std::pair<double,double>> in ascending order by second
void mysort::sort_vec_double2_y_increasing( std::vector<std::pair<double,double>> &vec_double2 )
{
  std::sort(vec_double2.begin(),vec_double2.end(),mysort::vec_double2_sort_y_increasing());
};
