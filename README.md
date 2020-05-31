

# Running the model

## Compiling and installing the NEST module

From the repository folder:

```
export NEST_INSTALL_DIR=/path/to/nest/
cd nest_modules
mkdir mb
cd mb
cmake -Dwith-nest=${NEST_INSTALL_DIR}/bin/nest-config ../IzhikevichCond
make
make install
```

If everything went right you can load the NEST module in PyNEST by doing `nest.InstallModule('izhikevich_cond_module')`.



## Running the model with normal behavior





- Parámetros, material suplementario y código
- el codigo fuente y el "experimental setup"
- vendría bien un README diciendo como descargarte el codigo, instalar los programas necesarios y como ejecutar el código
- ademas, este README debería tener una referencia al artículo publicado (eso se modificara una vez publica el articulo) y una pequeña descripción del modelo (capas del modelo, numero de neuronas en cada capa, tipos de sinapsis)
