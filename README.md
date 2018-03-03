# ProjetCrypto


Installation de Sodium :

Une fois dans le dossier libsodium-stable

./configure
make && make check
sudo make install

Compilation :

gcc -o elgamal elGamal.c -lgmp -lsodium


