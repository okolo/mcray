# Installing external libraries

<pre><code>mkdir SOPHIA
mkdir nr</code></pre>

- download SOPHIA event generator sources from https://elsevier.digitalcommonsdata.com/datasets/pkx5j87mgn/1 and extract in **SOPHIA** folder:
<pre><code>wget "https://prod-dcd-datasets-cache-zipfiles.s3.eu-west-1.amazonaws.com/pkx5j87mgn-1.zip"
unzip pkx5j87mgn-1.zip
cd SOPHIA/
tar -xzvf '../Monte Carlo simulations of photohadronic processes in astrophysics/adlb_v1_0.tar.gz'
cd ../
rm -r 'Monte Carlo simulations of photohadronic processes in astrophysics'
rm pkx5j87mgn-1.zip
</code></pre>
- download Numerical Recipes 3rd ed. C++ source code from http://numerical.recipes/com/storefront.html and extract to **nr** folder. TODO: fix this link
