$$
\frac{dX(t)}{dt} = c3 + c1 X\left( t \right) - c2 X\left( t \right)
$$

Physiological model:
$$
\begin{align}
\frac{ds(t)}{dt} =& \mathrm{tvr\_b}\left( t \right) - \mathrm{dec\_s}\left( t \right) \\
\mathrm{db}\left( t \right) =& Y \mathrm{u1}\left( t \right) - \mathrm{tvr\_b}\left( t \right) \\
\frac{db(t)}{dt} =& \mathrm{db}\left( t \right) \\
\frac{dcr_{tot(t)}}{dt} =& \mathrm{r\_tot}\left( t \right) \\
\frac{dr(t)}{dt} =& Y \mathrm{u1}\left( t \right) \left( b\left( t \right) \right)^{-1} \left(  - r\left( t \right) + s\left( t \right) \left( kmr + s\left( t \right) \right)^{-1} \right) \\
\mathrm{u1}\left( t \right) =& ks b\left( t \right) r\left( t \right) s\left( t \right) \left( km + s\left( t \right) \right)^{-1} \\
\mathrm{u2}\left( t \right) =& ms b\left( t \right) s\left( t \right) \left( km + s\left( t \right) \right)^{-1} \left( 1 - r\left( t \right) \right) \\
\mathrm{dec\_s}\left( t \right) =& \mathrm{u1}\left( t \right) + \mathrm{u2}\left( t \right) \\
\mathrm{tvr\_b}\left( t \right) =& kd b\left( t \right) \\
\mathrm{r\_gr}\left( t \right) =& \mathrm{u1}\left( t \right) \left( 1 - Y \right) \\
\mathrm{r\_m}\left( t \right) =& \mathrm{u2}\left( t \right) \\
\mathrm{r\_tot}\left( t \right) =& \mathrm{r\_gr}\left( t \right) + \mathrm{r\_m}\left( t \right) \\
q\left( t \right) =& \mathrm{u1}\left( t \right) \left( HB Y - HS \right) - HS \mathrm{u2}\left( t \right)
\end{align}
$$


