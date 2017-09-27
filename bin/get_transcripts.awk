#!/usr/bin/gawk
BEGIN {
    filefirst[array] = 1;
}
{
    match(FILENAME, /(.+)\/.+/, name);
    if(filefirst[FILENAME] != 0) {
        print $col >> outdir  name[1] "_" type ".tbl"
    }
    filefirst[FILENAME] = 1;
}
