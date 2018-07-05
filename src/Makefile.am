# ---------------------------------------------------------------------------
#
#   Kiri's Cropping Tool
#
#   Copyright (c) 2018 Kiri Jolly
#     http://expiredpopsicle.com
#     expiredpopsicle@gmail.com
#
# ---------------------------------------------------------------------------
#
#   This software is provided 'as-is', without any express or implied
#   warranty. In no event will the authors be held liable for any
#   damages arising from the use of this software.
#
#   Permission is granted to anyone to use this software for any
#   purpose, including commercial applications, and to alter it and
#   redistribute it freely, subject to the following restrictions:
#
#   1. The origin of this software must not be misrepresented; you must
#      not claim that you wrote the original software. If you use this
#      software in a product, an acknowledgment in the product
#      documentation would be appreciated but is not required.
#
#   2. Altered source versions must be plainly marked as such, and must
#      not be misrepresented as being the original software.
#
#   3. This notice may not be removed or altered from any source
#      distribution.
#
#
# -------------------------- END HEADER -------------------------------------

bin_PROGRAMS = kcrop

kcrop_SOURCES = \
	main.cpp \
	pixelimage_png.cpp pixelimage_png.h \
	pixelimage_jpeg.cpp pixelimage_jpeg.h

kcrop_CXXFLAGS = \
	$(libpng_CFLAGS) \
	$(libjpeg_CFLAGS) \
	$(liblilyengine_CFLAGS) \
	-Wall \
	--std=c++11 -O3 -g

kcrop_LDFLAGS = \
	$(libpng_LIBS) \
	$(libjpeg_LIBS) \
	$(liblilyengine_LIBS) \
	-Wall \
	 --std=c++11 -O3 -g



