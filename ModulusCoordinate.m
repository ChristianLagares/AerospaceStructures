function [E, a] = ModulusCoordinate(y, z, Modulus, alpha)
    % ModulusCoordinate: f64(f64,f64)
    % Returns the modulus of elasticity given the coordinate. The function
    % will generate results based on the conditions imposed for the INME
    % 4717 Final Test. 
    [b1, ~, b3, b4, b5, b6, c] = LinearNACA0012();
    if z >= 0 && z <= b1
        if y > 0
            E = Modulus(1);
            a = alpha(1);
        else
            E = Modulus(2);
            a = alpha(2);
        end
    elseif z > b1 && z <= (0.402+0.171)*c && abs(y) >= b6
        if y > 0
            E = Modulus(5);
            a = alpha(5);
        else
            E = Modulus(4);
            a = alpha(4);
        end
    elseif z > b1 && z <= (0.402+0.171)*c && abs(y) <= b6
        E = Modulus(6);
        a = alpha(6);
    elseif z < 0
        E = Modulus(3);
        a = alpha(3);
    elseif y > b6 | z < -0.128*c | z > (0.402+0.172)*c 
        error('Coordinates provided are outside the wingbox.')
    else
        error('Please check coordinates.')
    end
end