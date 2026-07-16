// tuple_int.cpp
#include "ns_myapp.hpp"
#include "ns_tuple_int.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"

//#########################################################################
//#########################################################################
// tuple_int functions
//#########################################################################
//#########################################################################
// equal for std::tuple<int,int>
bool tuple_int::equal( const std::tuple<int,int> &tp0, const std::tuple<int,int> &tp1)
{
  const int& i0_tp0 = std::get<0>(tp0);
  const int& i0_tp1 = std::get<0>(tp1);
  if( i0_tp0 != i0_tp1 ) return false;

  const int& i1_tp0 = std::get<1>(tp0);
  const int& i1_tp1 = std::get<1>(tp1);
  if( i1_tp0 != i1_tp1 ) return false;
  
  return true;
};

// output functions for UmpIntInt3
void tuple_int::out_ump_int_int2(FILE* fout, const UmpIntInt2& ump_int_int2)
{
  for (const auto& [key, val] : ump_int_int2) {
    const auto& [int0, int1] = val;
    fprintf(fout, "%d %d %d\n", key, int0, int1);
  }
}

void tuple_int::out_ump_int_int2(
  const std::filesystem::path& pathout, const UmpIntInt2 &ump_int_int2 )
{
  FILE *fout = myapp::get_fout(pathout);
  out_ump_int_int2(fout,ump_int_int2);
  myapp::close(fout,pathout);
}

void tuple_int::out_ump_int2_int( FILE *fout, const UmpInt2Int &ump_int2_int )
{
  for(const auto&[p_i0i1, i2]: ump_int2_int ){
    const auto & [i0, i1] = p_i0i1;
    fprintf(fout,"%d %d %d\n",i0,i1,i2);
  }
}

void tuple_int::out_ump_int2_int( const std::filesystem::path& pathout, const UmpInt2Int &ump_int2_int )
{
  FILE *fout = myapp::get_fout(pathout);
  out_ump_int2_int(fout,ump_int2_int);
  myapp::close(fout,pathout);
}

// get key from value
int tuple_int::get_key_from_value( const UmpIntInt2 &ump, const std::tuple<int,int> &value )
{
  auto it = std::find_if( ump.begin(), ump.end()
    , [&value]( const std::pair< int, std::tuple<int,int> > &p)
    {
      return equal(value,p.second);
    }
  );
  if( it == ump.end() ){
    THROW_ERROR("tuple_int::get_key_from_value: key was not found with value = {}", to_string(value));
  }
  return it->first;
}

// Get int2-tuple value from int key
// Returns tp_int2_notfound (-1,-1) if key not found
std::tuple<int,int> tuple_int::get_value_from_key( const UmpIntInt2 &ump, const int &key )
{
  auto found = ump.find(key);
  if(found == ump.end()) return tp_int2_notfound;
  return found->second;
}

// Get int value from int2-tuple key
// Returns int_notfound (-1) if key not found
int tuple_int::get_value_from_key( const UmpInt2Int &ump, const std::tuple<int,int> &key )
{
  auto found = ump.find(key);
  if(found == ump.end()) return int_notfound;
  return found->second;
}

// Get int value from two ints (constructs int2-tuple key)
int tuple_int::get_value_from_key( const UmpInt2Int &ump, const int i0, const int i1 )
{
  const std::tuple<int,int> key = std::make_tuple(i0,i1);
  return get_value_from_key(ump,key);
}


std::string tuple_int::to_string( const std::tuple<int,int> &tp )
{
  std::ostringstream oss;
  oss << std::get<0>(tp) << " " << std::get<1>(tp);
  return oss.str();
}

std::string tuple_int::to_string( const int i0, const int i1 )
{
  std::ostringstream oss;
  oss << i0 << " " << i1;
  return oss.str();
}

std::tuple<int,int> tuple_int::string_to_int2( const std::string &str )
{
  std::vector<std::string> vec_str = myapp::split(str,char_delim);
  if( vec_str.size()!=2 ) THROW_ERROR("vec_str.size()!=2");
  const int i0 = std::stoi( vec_str.at(0) );
  const int i1 = std::stoi( vec_str.at(1) );
  return std::make_tuple( i0, i1 );
}


//#########################################################################
//#########################################################################
// std::tuple<int,int,int> functions
//#########################################################################
//#########################################################################
// equal for std::tuple<int,int,int>
bool tuple_int::equal( const std::tuple<int,int,int> &tp0, const std::tuple<int,int,int> &tp1 )
{
  const int& i0_tp0 = std::get<0>(tp0);
  const int& i0_tp1 = std::get<0>(tp1);
  if( i0_tp0 != i0_tp1 ) return false;

  const int& i1_tp0 = std::get<1>(tp0);
  const int& i1_tp1 = std::get<1>(tp1);
  if( i1_tp0 != i1_tp1 ) return false;
  
  const int& i2_tp0 = std::get<2>(tp0);
  const int& i2_tp1 = std::get<2>(tp1);
  if( i2_tp0 != i2_tp1 ) return false;

  return true;
};

// output functions for UmpIntInt3
void tuple_int::out_ump_int_int3(FILE* fout, const UmpIntInt3& ump_int_int3)
{
  for (const auto& [key, val] : ump_int_int3) {
    const auto& [i0, i1, i2] = val;
    fprintf(fout, "%d %d %d %d\n", key, i0, i1, i2);
  }
}

void tuple_int::out_ump_int_int3( const std::filesystem::path& pathout, const UmpIntInt3 &ump_int_int3 )
{
  FILE *fout = myapp::get_fout(pathout);
  out_ump_int_int3(fout,ump_int_int3);
  myapp::close(fout,pathout);
}

void tuple_int::out_ump_int3_int( FILE *fout, const UmpInt3Int &ump_int3_int )
{
  for(const auto& elm: ump_int3_int ){
    const std::tuple<int,int,int> tp = elm.first;
    const auto& [i0, i1, i2] = tp;
    const int uni_id = elm.second;
    fprintf(fout,"%d %d %d %d\n",i0,i1,i2,uni_id);
  }
}

void tuple_int::out_ump_int3_int( const std::filesystem::path& pathout, const UmpInt3Int &ump_int3_int )
{
  FILE *fout = myapp::get_fout(pathout);
  out_ump_int3_int(fout,ump_int3_int);
  myapp::close(fout,pathout);
}

// get key from value
int tuple_int::get_key_from_value( 
    const UmpIntInt3 &ump
  , const std::tuple<int,int,int> &value )
{
  auto it = std::find_if( ump.begin(), ump.end()
    , [&value]( const std::pair< int, std::tuple<int,int,int> > &p)
    {
      return equal(value,p.second);
    }
  );
  if( it == ump.end() ){
    THROW_ERROR("tuple_int::get_key_from_value: key was not found with value = {}", to_string(value));
  }
  return it->first;
}

// Get int3-tuple value from int key
// Returns tp_int3_notfound (-1,-1,-1) if key not found
std::tuple<int,int,int> tuple_int::get_value_from_key( const UmpIntInt3 &ump, const int &key )
{
  auto found = ump.find(key);
  if(found == ump.end()) return tp_int3_notfound;
  return found->second;
}

// Get int value from int3-tuple key
// Returns int_notfound (-1) if key not found
int tuple_int::get_value_from_key( const UmpInt3Int &ump, const std::tuple<int,int,int> &key )
{
  auto found = ump.find(key);
  if(found == ump.end()) return int_notfound;
  return found->second;
}

// Get int value from three ints (constructs int3-tuple key)
int tuple_int::get_value_from_key( const UmpInt3Int &ump, const int i0, const int i1, const int i2 )
{
  const std::tuple<int,int,int> key = std::make_tuple(i0,i1,i2);
  return get_value_from_key(ump,key);
}

std::string tuple_int::to_string( const std::tuple<int,int,int> &tp )
{
  std::ostringstream oss;
  oss << std::get<0>(tp) << " " << std::get<1>(tp) << " " << std::get<2>(tp);
  return oss.str();
}

std::string tuple_int::to_string( const int i0, const int i1, const int i2 )
{
  std::ostringstream oss;
  oss << i0 << " " << i1 << " " << i2;
  return oss.str();
}

std::tuple<int,int,int> tuple_int::string_to_int3( const std::string &str )
{
  std::vector<std::string> vec_str = myapp::split(str,char_delim);
  if( vec_str.size()!=3 ) THROW_ERROR("vec_str.size()!=3");
  const int i0 = std::stoi( vec_str.at(0) );
  const int i1 = std::stoi( vec_str.at(1) );
  const int i2 = std::stoi( vec_str.at(2) );
  return std::make_tuple( i0, i1, i2 );
}

// @brief Extracts the set of keys from a std::map<int,int>
std::set<int> tuple_int::get_set_key( const std::map<int,int> &map )
{
  std::set<int> set_key;
  for(const auto& [key, val] : map){
    set_key.insert(key);
  }
  return set_key;
}

// @brief Extracts the set of values from a std::map<int,int>
std::set<int> tuple_int::get_set_value( const std::map<int,int> &map )
{
  std::set<int> set_value;
  for(const auto& [key, val] : map){
    set_value.insert(val);
  }
  return set_value;
}

// @brief Extracts the vector of keys from a std::map<int,int>
std::vector<int> tuple_int::get_vec_key( const std::map<int,int> &map )
{
  std::vector<int> vec_key;
  for(const auto& [key, val] : map){
    vec_key.push_back(key);
  }
  return vec_key;
}

// @brief Extracts the vector of values from a std::map<int,int>
std::vector<int> tuple_int::get_vec_value( const std::map<int,int> &map )
{
  std::vector<int> vec_value;
  for(const auto& [key, val] : map){
    vec_value.push_back(val);
  }
  return vec_value;
}