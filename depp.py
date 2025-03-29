import deepl

# Ваш API-ключ (замените на свой)
AUTH_KEY = "379fa3bc-7867-46d9-bace-750f6fd673b3:fx"

# Создание клиента DeepL
translator = deepl.Translator(AUTH_KEY)

# Путь к загруженному файлу
input_pdf = "/mnt/data/Факторизая.pdf"
output_pdf = "/mnt/data/Факторизая_перевод.pdf"

# Перевод документа
translator.translate_document_from_filepath(
    input_pdf, output_pdf, target_lang="RU"
)

print("Перевод завершен! Файл сохранен:", output_pdf)
