## Скрипты

Для установки необходимых библиотек для скриптов можно использовать команду
```
    pip install -r requirements.txt
```

* [generation.py](generation.py) - для генерации функции распределения тепла в начальный момент на основе [setting.ini](../initial/setting.ini)

* [plotX.py](plotX.py) - для построения сечения вдоль плоскости X графиков распределения тепла в начальный и конечный *(результат работы программы)* момент времени. 

      Параметры входящей строки:
        1. Путь до файла результата относительно папки result

*пример ввода:*
```
    plotX.py Kirill/euler.txt
```

* [plot4D.py](plot4D.py) - для построения четырёхмерных графиков распределения тепла в начальный и конечный *(результат работы программы)* момент времени. 

      Параметры входящей строки:
        1. Путь до файла настроек относительно папки initial
        2. Путь до файла функции в начальный момент времени относительно папки initial
        3. Путь до файла функции в конечный момент времени относительно папки result

*пример ввода:*
```
   plot4D.py setting.ini function.txt Kirill/result.txt
```
* [fault.py](fault.py) - для получение абсолютной и относильной погрешностей по результатам двух файлов. Построение графиков погрешностей в каждой точке по X.

      Параметры входящей строки:
        1. Путь до первого файла результата относительно папки result
        2. Путь до второго файла результата относительно папки result
        
 *пример ввода:*
```
    fault.py Kirill/euler.txt Kirill/implicit.txt
```
