#!/bin/bash

#####
#	PARSE COMMAND LINE
#####
function usage() {
    echo "Usage: $0 EXECUTABLE" 1>&2
    exit 1
}

# A program must be defined as the first argument
EXECUTABLE="$1"
if [ "$EXECUTABLE" == "" ]
then
    echo "$0: no executable defined" 1>&2
    usage
fi


#####
#	SET PATH TO EXECUTABLE
#####
# Check if path (full or relative) is set
DIRNAME=$(dirname $EXECUTABLE)
BASENAME=$(basename $EXECUTABLE)


# If the path to the executable is set, use it
if [[ "$EXECUTABLE" == */* ]]
then
    IS_EXECUTABLE_THERE=$(which $EXECUTABLE)
    if [ "$IS_EXECUTABLE_THERE" == "" ]
    then
        # Do not look for other versions of the executable if path is set explicitly
        # This will help avoid confusion with other versions lying around
        echo "$0: did not find \"$BASENAME\" in user-defined path: $DIRNAME: exiting" 1>&2
        echo > hdefexec.tmp
        exit 1
    else
        echo "$0: found \"$EXECUTABLE\""
        echo $EXECUTABLE > hdefexec.tmp
    fi
    exit
fi


# If path is not set, assume that this script is in /path/to/hdef/test/ (the test/ directory of the Hdef source)
# Look for the Hdef executable in the following places (in this order):
#     ../build      If Hdef is built in /path/to/hdef/build/, this will contain the most recently compiled executables
#     $PATH         Look for the codes found in the user's environment


# Check if the executable is in the relative directory ../build/
IS_EXECUTABLE_THERE=$(which $(dirname $0)/../build/$BASENAME | xargs dirname)
if [ "$IS_EXECUTABLE_THERE" != "" ]
then
    DIRNAME=$IS_EXECUTABLE_THERE
    echo "$0: found \"$BASENAME\" in $DIRNAME"
    echo $DIRNAME/$BASENAME > hdefexec.tmp
    exit
else
    echo "$0: could not find \"$BASENAME\" in ../build/" 1>&2
fi


# Check if the executable is in the user's PATH
IS_EXECUTABLE_THERE=$(which $BASENAME | xargs dirname)
if [ "$IS_EXECUTABLE_THERE" != "" ]
then
    DIRNAME=$IS_EXECUTABLE_THERE
    echo "$0: found \"$BASENAME\" in user's PATH folder: $DIRNAME"
    echo $DIRNAME/$BASENAME > hdefexec.tmp
    exit
else
    echo "$0: could not find \"$BASENAME\" in user's PATH" 1>&2
fi


# Hdef executable not found
echo "$0: unable to find Hdef executable \"$BASENAME\"; exiting" 1>&2
echo > hdefexec.tmp
exit 1
