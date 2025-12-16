Чтобы использовать в Юпитере нужно вот таким образом добавить в path 

```
import sys
from pathlib import Path

# Add your package directory to Python path
package_path = Path("/home/fedorov_de/MGX/viral-snake/")  # parent of your_package folder
sys.path.insert(0, str(package_path))
```

А потом как обычно 

```
from viral_snake.mgi_import import  validate_and_clean_ids, normalize_barcode_ids
from viral_snake.mgi_import import  rename_lane_columns, extract_lane_columns
from viral_snake.mgi_import import  create_symlinks_multi_lane
```