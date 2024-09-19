#!/bin/bash
cd solver/src
make $1
cd ../..

cd physics/src
make $1
cd ../..

cd amr/src
make $1
cd ../..

