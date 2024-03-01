def calc():
    history = []
    while True:
        print('asdf')
        x, y = (yield)
        if x == "h":
            print(history)
            continue
        result = x + y
        print(result)
        history.append(result)
        x,y = (yield)
        print('tot')

def ca():
    while True:
        print('fdsa')
        x, y = (yield)


c = calc()
b = ca()
print(type(c)) # <type 'generator'>

c.__next__() # Необходимая инициация. Можно написать c.send(None)
c.send((1,2)) # Выведет 3
b.__next__()
c.send((100, 30)) # Выведет 130
c.send((666, 0)) # Выведет 666
b.send((1,1))
c.send(('h',0)) # Выведет [3, 130, 666]
c.close() # Закрывем генератор
