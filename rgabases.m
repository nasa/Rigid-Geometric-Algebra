function [e0,e1,e2,e3,e4,e23,e31,e12,e43,e42,e41,e321,e412,e431,e423,e1234] = rgabases
%RGABASES Create full set of basis elements for RGA
e0 = rga('e0');
e1 = rga('e1'); e2 = rga('e2'); e3 = rga('e3'); e4 = rga('e4');
e23 = rga('e23'); e31 = rga('e31'); e12 = rga('e12');
e43 = rga('e43'); e42 = rga('e42'); e41 = rga('e41');
e321 = rga('e321'); e412 = rga('e412'); e431 = rga('e431'); e423 = rga('e423');
e1234 = rga('e1234');
end