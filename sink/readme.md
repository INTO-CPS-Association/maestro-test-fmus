# Dymola sources have been patched with:

```bash
#!/bin/bash
find . -name fmiCoSimFunctions_int.c -exec sed -i.bak s/"\"%s: derivative order %d is not supported\", order\[i\], label);"/"\"%s: derivative order %d is not supported\", label, order\[i\]);"/ {} \;
```

# Model description patched with 


```xml
<ModelStructure/>
```

