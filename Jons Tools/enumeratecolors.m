function [color] = enumeratecolors(index)
    switch index
        case 1
            color = 'black';
        case 2
            color = 'blue';
        case 3
            color = 'red';
        case 4
            color = 'green';
        case 5
            color = 'magenta';
    end
end