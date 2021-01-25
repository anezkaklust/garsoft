#!/bin/bash

sed -i 's+http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd+GDMLSchema/gdml.xsd+g' $(find . -maxdepth 1 -type f -name \*.gdml)

sed -i 's+http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd+../GDMLSchema/gdml.xsd+g' $(find . -maxdepth 1 -not -path . -type d -print0 | xargs -0 -I "{}"  find "{}" -name \*.gdml -type f)