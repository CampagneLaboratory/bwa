#!/bin/sh -e

test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.
#autoreconf --force -I /opt/local/share/aclocal --install --verbose "$srcdir"
# Need to set the following variable based on OS. MacOS needs the -I option to aclocal
OS=`uname -s`
if [ $OS = "Darwin" ]; then
  INCLUDE_OPTION="-I /opt/local/share/aclocal"
 else
  INCLUDE_OPTION=""
fi

rm -f config.cache
rm -fr autom4te.cache

autoheader
autoreconf -i
aclocal -I m4 $INCLUDE_OPTION
automake
autoconf
