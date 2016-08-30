aclocal -I m4
autoheader -f
automake -f --add-missing
autoconf -f
