От данной функции я хочу :
1) Подавать и распарсить файл minimap2.log. За это отвечает minimap2_browser.py load
2) Команда minimap2_browser.py stats должна выводить статистику, сколько каждому таксону соответствует строк в таблице, индексы и регионы
3) Команда python minimap2_browser.py interactive должна запрашивать имя вируса и отриосвывать выравнивание на контиги

Можно запускать так : 
1) Команда для подгрузки датасета
python3 alignment_plot.py load "/mnt/mgx/DATASETS/INTERNAL/VIROME/RUN3/co_assembly/megahit/*/contigs_formatted_minlen_*/minimap2/*/minimap2.log"
2) Статистика 
python3 alignment_plot.py stats
3) Исследовать только определенный таксон
python3 alignment_plot.py search "Norwavirus_beijiense"
4) Интерактивная часть сильно страдает :
Нужно чтобы сохранялся файл как название таксона_регион (то есть колонка prefix). Сейчас нормально сохраняется только при указании четкого --output
python3 alignment_plot.py show 0 --output norwavirus_plot.html