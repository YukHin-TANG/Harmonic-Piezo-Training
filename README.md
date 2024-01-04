# Harmonic-Piezo-Training

This project is created to train the dynamic model of piezo-material.

If we want to get the response of piezo-material under $$ v(t) = V_1 \cos(\omega t) + V_2 \sin(\omega t)$$, we can use the following equation to calculate the response.

=> 
$$ v(t) = V_1 \cdot \frac{e^{j \omega t}+e^{-j \omega t}}{2} + V_2 \cdot \frac{e^{j \omega t}-e^{-j \omega t}}{2j}= \frac{V_1-jV_2}{2} \cdot e^{j \omega t} + \frac{V_1+jV_2}{2} \cdot e^{-j \omega t}$$