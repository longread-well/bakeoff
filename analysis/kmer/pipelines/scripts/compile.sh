#!/bin/bash
/apps/well/gcc/7.2.0/bin/g++ --static -I /apps/well/boost/1.57.0/include/ -L /apps/well/boost/1.57.0/lib/ -lboost_iostreams -o make_mhap_kmer_file -lz /users/kwiatkowski/gav//Projects/Software/qctool/3rd_party/boost_1_55_0/libs/iostreams/src/zlib.cpp /users/kwiatkowski/gav//Projects/Software/qctool/3rd_party/boost_1_55_0/libs/iostreams/src/gzip.cpp /users/kwiatkowski/gav//Projects/Software/qctool/3rd_party/boost_1_55_0/libs/iostreams/src/file_descriptor.cpp -lz make_mhap_kmer_file.cpp