#### Uso de los programas:
1. Hacer el build de Mesher del siguiente modo:  
  1.1 en una carpeta correr: cmake [path a carpeta src]  
  1.2 en la misma carpeta: make
2. Correr Mesher para generar la malla de volumen con el comando:  
  ./mesher_roi -o [malla.off] -a [nivel de refinacion] -m
3. Obtener la superficie de la malla con Viewer  
  3.1 se debe compilar de la misma manera que Mesher  
  3.2 correr Viewer con el comando:
   ./viewer -m [malla.m3d]
Esto genera un archivo de GeomView que está en formato .off que contiene la superficie de la malla generada
4. Correr Mesher con el siguiente comando:
   ./mesher_roi -o [malla original.off] -o [malla generada.off] -a [nivel de refinacion del paso 2]


En este último paso se mostrará por pantalla los siguientes datos:  
- Distancia de Hausdorff desde malla original a generada
- Error de representación promedio
- Error de representación máximo
- Distancia de Hausdorff desde malla generada a original
- cantidad total de octantes de la superficie de la malla generada
- cantidad de octantes con un punto sobre el umbral de error

#### Consideraciones:  
- El umbral de error está hardcodeado porque no logré hacer que el programa lo considerara como argumento. En el código tiene un valor del 1% y para cambiarlo está en `main.cpp` en la línea 298. 
- La distancia de Hausdorff debería tener valores no superiores a 1, lo que ocurre la mayoría de las veces pero hay casos en que al correr el programa se obtiene un valor muy alto.
En ese caso se debería correr el programa nuevamente.
- Como mejora futura puede simplificar el uso de los programas implementando la obtención de la superficie dentro de Mesher, y así dejar de utilizar Viewer, y solo tener que correr Mesher una única vez.

#### Ejemplos  
se asume que el programa se corre desde un directorio build

./mesher_roi -o ../meshes/fertility_cgal.off -o ../fertility-off/fert2.off -a 2  
input->generated distance: -4.50289e-07  
average error: 0.0134537  
max error: 0.0422221  
generated->input distance: 0.842602  
total surface octants: 64  
octants over threshold: 58  

./mesher_roi -o ../meshes/fertility_cgal.off -o ../fertility-off/fert3.off -a 3  
input->generated distance: 1.08005e-06  
average error: 0.0068098  
max error: 0.0221073  
generated->input distance: 0.0441182  
total surface octants: 487  
octants over threshold: 210  

./mesher_roi -o ../meshes/protesis.off -o ../protesis-off/prot3.off -a 3  
input->generated distance: -3.27307e-07  
average error: 0.00411394  
max error: 0.0140368  
generated->input distance: 0.0259612  
total surface octants: 503  
octants over threshold: 55  
