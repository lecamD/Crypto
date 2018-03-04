# ProjetCrypto

Installation de GMP : 

https://gmplib.org/

Une fois dans le dossier gmp

./configure
make && make check
sudo make install

gcc -o elgamal elGamal.c -lgmp

Installation de Sodium :

https://download.libsodium.org/libsodium/releases/

Version du 2 mars 2018

sudo apt-get install libsodium-dev

Une fois dans le dossier libsodium-stable

./configure
make && make check
sudo make install

Compilation :

gcc -o elgamal elGamal.c -lgmp -lsodium


