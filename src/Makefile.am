
bin_PROGRAMS = intercrop

intercrop_SOURCES = \
	main.cpp \
	pixelimage_png.cpp pixelimage_png.h \
	pixelimage_jpeg.cpp pixelimage_jpeg.h

intercrop_CXXFLAGS = \
	$(libpng_CFLAGS) \
	$(libjpeg_CFLAGS) \
	$(liblilyengine_CFLAGS) \
	-Wall \
	--std=c++11

intercrop_LDFLAGS = \
	$(libpng_LIBS) \
	$(libjpeg_LIBS) \
	$(liblilyengine_LIBS) \
	-Wall \
	 --std=c++11



