#!/bin/bash
# $1 for type: remote_to_local or local_to_remote
# $2 local path
# $3 remote path
if [ $1 == 'local_to_remote' ]
then
    scp $2 sunchu@192.168.133.35:$3
elif [ $1 == 'remote_to_local' ]
then
    sudo scp -r sunchu@192.168.133.35:$3 $2   
else
    printf "************************************************************\n"
    printf "* Error: please provide remote_to_local or local_to_remote.*\n"
    printf "************************************************************\n"
    exit 1
fi