
%:
	g++ $@.cpp \
		-std=c++14 \
		-O3 \
		-DNDEBUG \
		-pedantic \
		-fopenmp \
		-lpthread \
		-lrt \
		-lz \
		-DSEQAN_HAS_ZLIB=1 \
		-DSEQAN_ENABLE_DEBUG=0 \
		-DSEQAN_ENABLE_TESTING=0 \
		-o $@

