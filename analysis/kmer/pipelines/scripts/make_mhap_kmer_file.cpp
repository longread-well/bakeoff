#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

std::vector< char > makeReverseComplementTable() {
	std::vector< char > result( 256 ) ;
	for( std::size_t i = 0; i < 256; ++i ) {
		result[i] = char(i) ;
	}
	result[ unsigned('A')] = 'T' ;
	result[ unsigned('a')] = 't' ;
	result[ unsigned('T')] = 'A' ;
	result[ unsigned('t')] = 'a' ;
	result[ unsigned('G')] = 'C' ;
	result[ unsigned('g')] = 'c' ;
	result[ unsigned('C')] = 'G' ;
	result[ unsigned('c')] = 'g' ;
	return result ;
}

void reverse_complement( std::string& kmer, std::vector< char > const& table ) {
	if( kmer.size() == 0 ) {
		return ;
	}
	std::size_t i = 0 ;
	std::size_t j = kmer.size() - 1 ;
	for( ; i < j; ++i,--j ) {
		char x = kmer[i] ;
		kmer[i] = table[ unsigned(kmer[j]) ] ;
		kmer[j] = table[ unsigned(x) ] ;
	}
	if( i == j ) {
		kmer[i] = table[ unsigned(kmer[i]) ] ;
	}
}

int main( int argc, char** argv ) {
	if( argc < 2 ) {
		std::cerr << "Supply an argument, the name of a file to load.\n" ;
		exit(-1) ;
	}

	long count = 0 ;
	long total = 0 ;	
	// first pass
	std::cerr << "Pass 1...\n" ;
	{
		std::ifstream file( argv[1], std::ios_base::in | std::ios_base::binary );
		boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
		in.push( boost::iostreams::gzip_decompressor() );
		in.push( file );
		std::istream inf(&in);

		std::string line ;
		while( std::getline( inf, line ) ) {
			std::size_t where = line.find( ' ' ) ;
			if( where == std::string::npos ) {
				std::cerr << "Input malformed on this line:\n"
					<< line << "\n"
					<< "...quitting.\n" ;
				exit(-1) ;
			}
			std::string kmerCount = line.substr( where + 1, line.size() ) ;
			total += std::stol( kmerCount ) ;
			++count ;
			if( count % 10000000 == 0 ) {
				std::cerr << "...read " << count << "...\n" ;
			}
		}
	}
	std::cerr << "Read " << count << " kmers with total coverage " << (double(total)/count) << "...\n" ;

	std::cerr << "Pass 2...\n" ;
	std::cout << count << "\t" << count << "\n" ;
	{
		std::vector< char > reverse_complement_table = makeReverseComplementTable() ;
		std::ifstream file( argv[1], std::ios_base::in | std::ios_base::binary );
		boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
		in.push( boost::iostreams::gzip_decompressor() );
		in.push( file );
		std::istream inf(&in);

		std::string line ;
		long outputCount = 0 ;
		while( std::getline( inf, line ) ) {
			std::size_t where = line.find( ' ' ) ;
			if( where == std::string::npos ) {
				std::cerr << "Input malformed on this line:\n"
					<< line << "\n"
					<< "...quitting.\n" ;
				exit(-1) ;
			}
			std::string kmer = line.substr( 0, where ) ;
			std::string kmerCount = line.substr( where + 1, line.size() ) ;
			std::cout << kmer << "\t" << std::fixed << std::setprecision(8) << double( std::stol( kmerCount ) ) / total  << "\n" ;
			reverse_complement( kmer, reverse_complement_table ) ;
			std::cout << kmer << "\t" << std::fixed << std::setprecision(8) << double( std::stol( kmerCount ) ) / total  << "\n" ;
			++outputCount ;
			if( outputCount % 10000000 == 0 ) {
				std::cerr << "...wrote " << outputCount << "...\n" ;
			}
		}
	}

	std::cerr << "Output " << count << " kmers with total coverage " << (double(total)/count) << "...\n" ;
	std::cerr << "Thank you for using make_mhap_kmer_file!\n" ;
}

