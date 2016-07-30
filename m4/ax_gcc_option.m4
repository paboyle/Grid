AC_DEFUN([AX_GCC_OPTION], [
AC_REQUIRE([AC_PROG_CC])

AC_MSG_CHECKING([if gcc accepts $1 option])

AS_IF([ test "x$GCC" = "xyes" ],[
AS_IF([ test -z "$3" ],[
ax_gcc_option_test="int main()
{
return 0;
}"
],[
ax_gcc_option_test="$3"
])

# Dump the test program to file
cat <<EOF > conftest.c
$ax_gcc_option_test
EOF

# Dump back the file to the log, useful for debugging purposes
AC_TRY_COMMAND(cat conftest.c 1>&AS_MESSAGE_LOG_FD)

AS_IF([ AC_TRY_COMMAND($CC $2 $1 -c conftest.c 1>&AS_MESSAGE_LOG_FD) ],[
AC_MSG_RESULT([yes])
$4
],[
AC_MSG_RESULT([no])
$5
])
],[
AC_MSG_RESULT([no gcc available])
])
])
