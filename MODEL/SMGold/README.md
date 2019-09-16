# Table of Contents

1.  [\(\phi^+ \rightarrow A_\mu \,\phi^-\)](#org3c23570)



<a id="org3c23570"></a>

# $`\phi^+ \rightarrow A_\mu \,\phi^-`$

![img](ltximg/eq75.png "EQ75")

This is loaded in the model file as [Model.C](Model.C):
``` c++
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_phiplus));
    m_v.back().AddParticle(Flavour(kf_phiplus).Bar());
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("GGV");
    m_v.back().cpl.push_back(-I*g1);
    m_v.back().order[1]=1;
```
The default ordering corresponds thus to [SSV_LC.C](SSV_LC.C) 
``` c++
    if (p_v->V()->id.back()==2)
```
which implies

-   $`\phi^+ \to 0 \, (a)`$
-   $`\phi^- \to 1 \, (b)`$
-   $`A_\mu  \to 2`$

The other two cases are, respectively:
``` c++
    if (p_v->V()->id.back()==1)
```
which implies

-   $`\phi^+ \to 1 \,(b)`$
-   $`\phi^- \to 2`$
-   $`A_\mu  \to 0 \,(a)`$

and
``` c++
    if (p_v->V()->id.back()==0)
```
which implies

-   $`\phi^+ \to 2`$
-   $`\phi^- \to 0 \,(a)`$
-   $`A_\mu  \to 1 \,(b)`$

In all cases, the resulting vertex is a scalar. Either because you
get the photon multiplied by its pol vector and a scalar or because
there are two scalars outgoing.

