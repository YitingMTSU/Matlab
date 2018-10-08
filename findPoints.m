function y = findPoints(f, a, b)
sym x;
x1 = a;
x2 = b;
eps = 1e-4;
for i = 1:3000
    x = x1;
    y1 = eval(f);
    x = x2;
    y2 = eval(f);
    if abs(y1)<eps
        y = x1;
    elseif abs(y2)<eps
        y = x2;
    else
        temp = (x1 + x2)/2;
        x = temp;
        y_temp = eval(f);
        if abs(y_temp)<eps
            y = temp;
        elseif y_temp * y1 > 0
            x1 = temp;
        else
            x2 = temp;
        end
    end
end
if i == 2000
    y = x1;
    i
end