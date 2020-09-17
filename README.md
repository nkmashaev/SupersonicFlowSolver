[a relative link](Description.md)

# Библиотека PSST для расчета сверхзвукового течения

## Описание
- - - - - - - - - - - - -

Программа разрабатывалась для запуска на кластере(Linux). Данный модуль расчитывает сверхзвуковое течения газа на неструктурированной сетке. Формат сетки соответствует формату хранения реализованном в коммерческом пакете ANSYS Fluent.  Для запуска программы также необходимо создать в текстовом редакторе следующие файлы: problem_manager, \*.bvp, \*.btq, \*.inp, \*.par, \*.prop. 

[Описание формата сетки](https://www.afs.enea.it/project/neptunius/docs/fluent/html/ug/node1464.htm#format-grid)

Примеры файлов со входными данными: 
- [problem_manager](./problem_manager);
- [wedge.inp](./wedge.inp); 
- [wedge.prop](./wedge.prop); 
- [wedge.par](./wedge.par); 
- [wedge.btq](./wedge.btq);
- [wedge.bvp](./wedge.bvp);
- [wedge.msh](./wedge.msh);
