## Downloading Conventional Atm Obs

from within this directory, rum `make` to compile  the following:

* `util/grabbufr` - NCEP utility to decode downloaded BUFR observation files
* `util/dec_prepbufr` - converts downloaded PREPBUFR files into a format readable by the LETKF

To download conventional atmospheric observations for a month, run:
```
./util/get_ncepobs.sh EMAIL PASSWD 20100101 20100131
```
where `EMAIL` and `PASSWD` are for your account at rda.ucar.edu