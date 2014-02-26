mkdir pack_bin
del /Q /S pack_bin\*.*
mkdir pack_bin\50
mkdir pack_bin\52
mkdir pack_bin\54
mkdir pack_bin\56
copy Release\ffans.exe pack_bin\50\
copy Release\ffans.exe pack_bin\52\
copy Release\ffans.exe pack_bin\54\
copy Release\ffans.exe pack_bin\56\
copy *param* pack_bin\50\
copy *param* pack_bin\52\
copy *param* pack_bin\54\
copy *param* pack_bin\56\