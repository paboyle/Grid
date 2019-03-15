#!/usr/bin/env bash

for m in `find Modules -name '*.hpp' -type f -print`; do
    echo "====== ${m}"
    CCFILE=`echo $m | sed -E s/\.hpp/.cc/g`
    NS=`echo $m | awk -F '/' '{print $2}'`
    NMOD=`grep -E 'MODULE_REGISTER_TMP.+<.+>.?' $m | wc -l`
    if [ ! -e ${CCFILE} ] && (( NMOD != 0 )); then
        echo "#include <Grid/Hadrons/${m}>" >> ${CCFILE}
        echo "" >> ${CCFILE}
        echo "using namespace Grid;" >> ${CCFILE}
        echo "using namespace Hadrons;" >> ${CCFILE}
        echo "using namespace ${NS};" >> ${CCFILE}
        echo "" >> ${CCFILE}
        for i in `grep -E 'MODULE_REGISTER_TMP.+<.+>.?' $m | sed -E 's/ +//g'`
        do
            TMPARG=`echo ${i} | grep -oE 'ARG\(.+>\)' | sed -E 's/^ARG\(//g' | sed -E 's/\)$//g'`
            SUB=`echo ${i} | sed -E 's/ARG\(.+>\)/@arg@/g' | sed -E 's/,/|/g'`
            SUB=`echo ${SUB} | sed -E 's/.+\(//g' | sed -E 's/\);//g'`
            SUB=`echo ${SUB} | sed -E "s/@arg@/${TMPARG}/g"`
            NAME=`echo ${SUB} | awk -F '|' '{print $1}'`
            TYPE=`echo ${SUB} | awk -F '|' '{print $2}'`
            echo "template class Grid::Hadrons::${NS}::${TYPE};" >> ${CCFILE}
        done
        echo "" >> ${CCFILE}
    fi
done