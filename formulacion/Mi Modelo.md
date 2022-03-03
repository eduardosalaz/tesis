## Modelo modificado para las MFI

### P-centro vs P-mediana:
P-centro minimiza la distancia máxima de viaje a su centro P más cercano. Por su parte, P-mediana minimiza el promedio de las distancias de los centros de demanda a los p centros.
Se usará una función de dispersión P-mediana.
revisar que tipo de función se usa para el riesgo
ignorar similitud
ignorar asignacion conjunta/disconjunta
branches de distintos tipos
penalizar riesgo
se puede modelar como una funcion para minimizar riesgo
y que la distancia se puede poner como una restriccion en lugar de la funcion objetivo
ir formulando de manera combinatoria el modelo

## Modelo combinatorio

Contamos con un conjunto $B$ en el que existen las unidades básicas de clientes (BUs) que hay que servir.
Adicionalmente existe un conjunto $S$ de $s$ sucursales, las cuales pueden pertenecer a uno de los 5 $k$ tipos de sucursales de acuerdo al parámetro $s_{ik}$ es igual a 1 si la sucursal $i$ es del tipo $k$ para todo $k \in 1\ldots5$, 0 de otro modo; tomando en cuenta que deben de ser repartidas de acuerdo a cotas inferiores y superiores de tipos de sucursales $l_k$ y $u_k$. 
Tenemos una matriz de distancias $d_{ij}$ para cada $i$ BU y una $j$ sucursal correspondiente.
Contamos con los siguientes conjuntos de variables de decisión binarias:
- $X_{ij}$, 1 si asignamos a la $i$ sucursal la $j$ BU, 0 si no, para toda $i \in S, j \in B$
- $Y_i$, 1 si usamos la sucursal $i$ como un centro territorial, 0 si no, para todo $i \in S$ 
Podemos hacer uso de distintas funciones estadísticas para determinar el riesgo asociado con ofrecer créditos a cada BU tomando en cuenta información histórica de sus ingresos. Varianza, desviación estándar, semivarianza o Valor en riesgo. Sea cual sea la función que usemos para determinar el riesgo y su balanceo, asociamos el parámetro $R_j$ para cada $j$ BU representando su propio riesgo.
Para poder balancear los terrenos en base al riesgo, asignamos un umbral de riesgo $\gamma_i$ para cada $i$ sucursal, tomando en cuenta que cada sucursal es el centro territorial.
Medimos 3 actividades $m$ para cada $j$ BU con el parámetro $v_j^m$, y asignamos una tolerancia territorial para cada medida de actividad $t^m$.

## Función objetivo
Minimizar 
$$ 
\sum_{i \in S} \sum_{j \in B} X_{ij}d_{ij} 
$$
## Restricciones
$$
\sum_{i \in S} X_{ij} = 1, \forall j \in B
$$
Asignar sólo a una sucursal $i$ la BU $j$

$$
X_{ij} \le Y_i, \forall i \in S, j \in B
$$
Sólo se pueden asignar a las sucursales determinadas como centros territoriales
$$
Y_i\mu_m^i(1-t^m) \le \sum_{j\in B}X_{ij}v_j^m \le Y_i\mu_m^i(1+t^m), \forall i \in S, m = 1\ldots3
$$
Medidas de actividad para cada territorio que estén dentro de un rango de tolerancia territorial.

$$
l_k \le \sum_{i \in S_k} Y_i \le u_k, k = 1 \ldots 5
$$
Asegurarse que se respeten los límites inferiores y superiores $l_k, u_k$ de los tipos de sucursales a emplearse de acuerdo a $k$. Definir conjuntos $s_k$ para cada sucursal que pertenezca a $k$, podemos usar cada conjunto $S_k$ en lugar de un parámetro binario. Hay que mantener una coherencia de los límites con respecto a la cardinalidad de cada $S_k$ ya sea ingresado por el usuario o calculado en base a ese tamaño.
Las sumatorias de las cotas deben de ser menores y mayores a todas las que tenemos de manera respectiva. El $u_k$ debe de ser menor que la cardinalidad de $S_k$, por su parte la sumatoria de $l_k$ debe de ser menor que la cardinalidad de $S$ y la sumatoria de las $u_k$ debe de ser mayor que la cardinalidad de $S$.
$$
\sum_{i \in S} Y_i = p
$$
El número de centros de territoriales debe concordar con el parámetro p

$$
\sum_{j \in B}X_{ij}R_j \le \alpha_i, \forall i \in S
$$
Para balancear el riesgo hay que escribir la sumatoria del parámetro de riesgo de cada BU por la $X_{ij}$ asignadas a la sucursal $i$ es la medida del riesgo del territorio $i$ acotar a una constante superior $\alpha_i$ para cada $i$ que también debe de tener sentido con respecto a los $R_j$ para la factibilidad del problema. Experimentar en base a los $R$.

Ir creando instancias (primero a mano, luego hago el generador) y su programación matemática correspondiente en CPLEX/Gurobi.
