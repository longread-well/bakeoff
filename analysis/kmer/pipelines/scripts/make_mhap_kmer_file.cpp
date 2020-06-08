#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

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
	std::cout << count << "\t" << 1024 << "\n" ;
	{
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
			++outputCount ;
			if( outputCount % 10000000 == 0 ) {
				std::cerr << "...wrote " << outputCount << "...\n" ;
			}
		}
	}

	std::cerr << "Output " << count << " kmers with total coverage " << (double(total)/count) << "...\n" ;
	std::cerr << "Thank you for using make_mhap_kmer_file!\n" ;
}

