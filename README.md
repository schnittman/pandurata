This is the version of Pandurata that was used in Schnittman, Krolik, and Noble (2013) [ApJ 769, 156]
It uses the IDL routine read_harm3d_data2 to process the output data from Harm3d found in /rawdata, re-writing it in a Pandurata-friendly form in /data

in /data, we need to convert the ascii data to binary via 
> gcc -o write_harm3d_bin write_harm3d_bin.c -lm
> ./write_harm3d_bin

compile Pandurata in the main directory via
> make -f makepan
> ./Pandurata
This will generate a bunch of data files in /data
