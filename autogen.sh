#!/bin/sh

if [ \! -e Lily-Engine-Utils ]; then
    git clone https://github.com/ExpiredPopsicle/Lily-Engine-Utils.git
fi

aclocal && \
automake --add-missing && \
autoconf
