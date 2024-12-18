#!/bin/bash

sudo mount -t cifs //192.168.2.50/Bioinformatics /mnt/nas -o credentials=~/.nas-credentials,vers=3.0,uid=1000,gid=1000
