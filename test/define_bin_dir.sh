#!/bin/bash

HOME_DIR=`echo $HOME`
cd ..
HDEF_DIR=`pwd`
cd test

grep "^ *BIN_DIR *=" ../Makefile |\
    tail -1 |\
    sed -e "s/.*=//" |\
    awk '{print $1}' |\
    awk -F/ '{
        if ($1==".") {
            if (NF==1) {
                printf("'$HDEF_DIR'")
            } else {
                printf("'$HDEF_DIR'/")
                for (i=2;i<=NF-1;i++) {
                    printf("%s/"),$i
                }
                printf("%s"),$NF
            }
        } else if ($1=="..") {
            if (NF==1) {
                printf("'$HDEF_DIR'/..")
            } else {
                printf("'$HDEF_DIR'/../")
                for (i=2;i<=NF-1;i++) {
                    printf("%s/"),$i
                }
                printf("%s"),$NF
            }
        } else if ($1=="~") {
            if (NF==1) {
                printf("'$HOME_DIR'")
            } else {
                printf("'$HOME_DIR'/")
                for (i=2;i<=NF-1;i++) {
                    printf("%s/"),$i
                }
                printf("%s"),$NF
            }
        } else {
            print $0
        }
    }'
