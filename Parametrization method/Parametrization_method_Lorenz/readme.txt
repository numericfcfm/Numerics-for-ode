Archivos de la carpeta Parametrization_method_Lorenz

Precaución. Antes de ejecutar cualquier archivo lea el articulo, El metodo de parametrizacion para variedades invariantes de
puntos de equilibrio de ecuaciones diferenciales ordinarias. Abstraction &
Application, 30. pp 64-81, (2020), que se encuentra en el repositorio.

1. El archivo parametrization_lorenz.py contiene un programa escrito en Python que calcula
la variedad estable (local) de dimensión dos asociada al punto de equilibrio (0,0,0) del sistema de Lorenz con los 
los parámetros "clasicos" donde surge el atractor de Lorenz. El programa viene precargado para calcular la variedad
con orden 20 y utilizando como dominio fundamental un disco de radio 16.  

2. El archivo global_manifold.m contiene código para Matlab donde se presenta un esquema básico para calcular
la variedad estable global del sistema de Lorenz. Para ejecutar el archivo necesita poner los archivos 
lorenzfield_normalized.m (campo vectorial con tiempo reparametrizado) y border.m (discretización de la frontera
de la variedad local del punto 1) en una misma carpeta. Para más referencias sobre este cálculo puede consultar 
la sección 2.4.2 de libro 

A. Haro, et.al. The Parametrization Method for Invariant Manifolds. From Rigorous Results to
Effective Computations. Applied Mathematical Sciences, Springer, 2016.