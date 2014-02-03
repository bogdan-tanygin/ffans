mkdir pack_bin
del /Q /S pack_bin\*.*
mkdir pack_bin\1
mkdir pack_bin\2
mkdir pack_bin\3
mkdir pack_bin\4
copy Release\ff-heating-exp-1.exe pack_bin\1\
copy Release\ff-heating-exp-1.exe pack_bin\2\
copy Release\ff-heating-exp-1.exe pack_bin\3\
copy Release\ff-heating-exp-1.exe pack_bin\4\
copy *param* pack_bin\1\
copy *param* pack_bin\2\
copy *param* pack_bin\3\
copy *param* pack_bin\4\